# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a 
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.2

###############################################################################
###########################   Imported modules   ##############################

import pandas as pd
import time
from numpy import arange, array


from Network_builder_public import create_pricesDF, Create_loads, heuristic_control

###############################################################################

def Load_data(P_L, start_day = 0, end_day = 364, dt = 0.25, H = 0):
    import csvreader    


    
    t_0 = start_day
    t_simulation = int((end_day - start_day)*24/dt)
    
    
    print('Loading data')
    CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
    CSVDataTamb = csvreader.read_data(csv='Tamb_15min.csv', address='')
    CSVDataRad = csvreader.read_data(csv='Radiation_1min.csv', address='')
    CSVDataPrices = csvreader.read_data(csv='DA_Prices_15min.csv', address='')
    CSVDataTsoil = csvreader.read_data(csv='Soil_dy_10cm.csv', address='')
    CSVDataPV.data2array()
    CSVDataTamb.data2array()
    CSVDataRad.data2array()
    CSVDataPrices.data2array()
    CSVDataTsoil.data2cols()
    P_PV_av = [i[0]/1000 for i in CSVDataPV.ar]
    T_amb = [i[0]+273 for i in CSVDataTamb.ar[t_0:t_0+t_simulation+H]]
    a = arange(0,len(CSVDataRad.ar),15)
    G = array([CSVDataRad.ar[i][0]/3600 for i in a[t_0:t_0+t_simulation+H]])
    Energy_price = [i[0]*dt for i in CSVDataPrices.ar[t_0:t_0+t_simulation+H]] # [energy_costs[-1]*0.25 for i in CSVDataPrices.ar]
    T_soil = CSVDataTsoil.col
    
    depth = 0.2
    dy = 0.1
    Height_TESS = 1.8
    T_soil = [i[int(depth/dy):int((depth + Height_TESS)/dy)] for i in T_soil]
    
    T_set_day = [273 + 20 - 3]*int((6-0)*4) + [273 + 20]*int((22-6)*4) + [273 + 20 - 3]*int((24-22)*4)
    T_set = array(T_set_day*(end_day+H-start_day))  
    
    df_Prices = create_pricesDF(Energy_price, t_simulation)   
    
    # P_L = get_base_load(node)
    
    P_PV_av = [int(sum(P_L)/sum(P_PV_av))*i for i in P_PV_av]    

    return([T_amb, T_set, T_soil, G, P_PV_av, df_Prices])


def get_base_load(house, start_day = 0, end_day = 364, dt = 0.25, H = 0):

    if house == 'studio':
        node = 46
    elif house == 'apartment':
        node = 5
    else:
        node = 71

    t_0 = start_day
    t_simulation = int((end_day - start_day)*24/dt)    
    
    DF_Network = pd.read_excel('Stedin_network_3.xlsx') # ['Gaia_network_2.xlsx', 'Stedin_network_2.xlsx', 'Stedin_network_3']
    Network_headers = ['Van.Nummer', 'Naar.Nummer', 'Lengte', 'Kabeltype', 'GM', 'House_Type']
    csv_name = 'S_n_Stedin_301.csv' # [S_n_Stedin_301_Case_2, 'S_n_CIGRE_18.csv', 'S_n_Stedin_301.csv']
        
    
    S_n = Create_loads(DF_Network, Network_headers, Load_data = True, csv_name = csv_name)    
    labels = [str(i+1) for i in range(len(S_n[0]))]
    df_P_L = pd.DataFrame(S_n[t_0:t_0+t_simulation+H], columns = labels)
        
    return df_P_L[str(node)].tolist()

def simulate_house(node, T_amb, T_set, T_soil, G, P_PV_av, P_L, df_Prices, House_type = 'apartment', enable_PV = False, enable_HP = False, enable_TESS = False, enable_BESS = False, start_day = 0, end_day = 364, dt = 0.25, P_ref = 10, DSO_operation = False, force_temperature = False, glazing_type = 'double', Cavity_Filling = False, Cavity_type = False, Wall_Cover = False, LRoof = 0.2, LCavity_wall = 0.05):

    if not enable_BESS:
        P_BESS_max = 0
    else:
        P_BESS_max = 10
    Capacity_BESS_BoL = 10
    enable_HP_2_TESS = enable_TESS

    t_simulation = int((end_day - start_day)*24/dt)
    
    # P_L = Load_data(node)
    
    # [T_amb, T_set, T_soil, G, P_PV_av, P_L, df_Prices] = Load_data(node)
    
    P_PV = [0 for i in range(t_simulation)]
    P_BESS = [0 for i in range(t_simulation)]
    P_HP = [0 for i in range(t_simulation)]
    P_HP_TESS = [0 for i in range(t_simulation)]
    P_Grid = [0 for i in range(t_simulation)]
    
    SOC_BESS = [0.5 for i in range(t_simulation)]
    C_BESS = [Capacity_BESS_BoL for i in range(t_simulation)]
    
    Qdot_D = [0 for i in range(t_simulation)]
    Qdot_HP = [0 for i in range(t_simulation)]
    Qdot_HP_TESS = [0 for i in range(t_simulation)]    
    Qdot_TESS = [0 for i in range(t_simulation)]
    Qdot_TESS_SD = [0 for i in range(t_simulation)]    
    Qdot_Boiler = [0 for i in range(t_simulation)]
    
    T_in = [20 + 273 for i in range(t_simulation)]
    T_TESS = [50 + 273 for i in range(t_simulation)]


    start_time = time.time()    # The timer is initializad.    
    for ts in range(1, t_simulation-1):
        # if (ts)%(24/dt) == 0:
            # print('Day', start_day/(24/dt) + int(ts/(24/dt)))        
        
        [P_PV[ts], P_BESS[ts], P_HP[ts], P_HP_TESS[ts], P_Grid[ts], SOC_BESS[ts], C_BESS[ts], Qdot_D[ts], Qdot_HP[ts], Qdot_HP_TESS[ts], T_TESS[ts], Qdot_TESS[ts], Qdot_TESS_SD[ts], Qdot_Boiler[ts], enable_HP_2_TESS_dummy, T_in[ts]] = heuristic_control(T_in[ts-1], T_set[ts], T_amb[ts], T_TESS[ts-1], 38+273, T_soil[ts], 0, P_L[ts], P_ref, DSO_operation, SOC_BESS[ts-1], df_Prices['Current'][ts], df_Prices['Median'][ts], df_Prices['Q25'][ts], df_Prices['Q75'][ts], C_BESS[ts-1], P_PV_av[ts], enable_HP = enable_HP, enable_TESS = enable_TESS, enable_HP_2_TESS = enable_HP_2_TESS, external_control = False, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, House_type = House_type, force_temperature = force_temperature, glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)
    end_time = time.time()    # The timer is finished
    
    # print('Computation time: ', (end_time - start_time)/60)
    
    return([P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in])

def Cost_windows(glazing_type, windows_area):
    if glazing_type == 'single':
        return 0
    
    elif glazing_type == 'double':
        return windows_area*700
        
    elif glazing_type == 'triple':
        return windows_area*1000

def Cost_roof(LRoof, roof_Area):

    if LRoof < 0.1:
        return 0
    elif LRoof <= 0.1:
        return roof_Area*35
    elif LRoof <= 0.15:
        return roof_Area*40
    elif LRoof <= 0.2:
        return roof_Area*50            
    else:
        return roof_Area*60
    
def Cost_Walls(walls_area, Cavity_Filling, Cavity_type, Wall_Cover):
    
    wall_cost = 0
    
    if Cavity_type == 'air':
        wall_cost += 0
    else:
        wall_cost += walls_area*40
        
    if Wall_Cover == False:
        wall_cost += 0
    else:
        wall_cost += walls_area*150
        
    return wall_cost
     

def Calculate_gas_consumption(Qdot_Boiler, CV_Gas = 12800, dt = 0.25):
    return dt*sum(Qdot_Boiler)/12800

def Calculate_thermal_demand(Qdot_Boiler, Qdot_HP, Qdot_TESS, dt = 0.25):
    # return dt*sum([-i for i in Qdot_D if i<0])/1e6
    return dt*(sum(Qdot_Boiler) + sum(Qdot_HP) + sum(Qdot_TESS))/1e6

def CAPEX(house_type, glazing_type, LRoof, Cavity_Filling, Cavity_type, Wall_Cover, enable_HP = False, PV_size = 0):
    from MCES_library import House_type_adjust
    
    [window_adjust, walls_adjust, roof_adjust, volume_adjust] = House_type_adjust(house_type)
    
    windows_area  = 8*window_adjust 
    roof_Area = 120.3*roof_adjust
    walls_area = 111.6*walls_adjust
    
    window_cost = Cost_windows(glazing_type, windows_area)   
    roof_cost = Cost_roof(LRoof, roof_Area)    
    walls_cost = Cost_Walls(walls_area, Cavity_Filling, Cavity_type, Wall_Cover)
    
    return window_cost + roof_cost + walls_cost + enable_HP*10000 + PV_size*1400

def OPEX(P_Grid, Qdot_Boiler, df_Prices, c_gas = 1.44, dt = 0.25):
    return dt*sum([p*c for p,c in zip(P_Grid, df_Prices['Current'].tolist())]) + Calculate_gas_consumption(Qdot_Boiler)*c_gas

def Energy_label_parameters(energy_label):
    
    if energy_label == 'A':
        glazing_type = 'triple'
        LRoof = 0.3
        LCavity_wall = 0.1
        Cavity_Filling = True
        Cavity_type = 'EPS'
        Wall_Cover = True
        
        enable_PV = True
        enable_HP = True
        enable_TESS = False
        enable_BESS = False

    elif energy_label == 'B':
        glazing_type = 'double'
        LRoof = 0.2
        LCavity_wall = 0.08
        Cavity_Filling = True
        Cavity_type = 'EPS'
        Wall_Cover = False 
        
        enable_PV = False
        enable_HP = False
        enable_TESS = False
        enable_BESS = False

    elif energy_label == 'C':
        glazing_type = 'double'
        LRoof = 0.15
        LCavity_wall = 0.08 
        Cavity_Filling = True
        Cavity_type =  'air'
        Wall_Cover = False 
        
        enable_PV = False
        enable_HP = False
        enable_TESS = False
        enable_BESS = False

    elif energy_label == 'D':
        glazing_type = 'double'
        LRoof = 0.15
        LCavity_wall = 0.05
        Cavity_Filling = True
        Cavity_type = 'air'
        Wall_Cover = False  
        
        enable_PV = False
        enable_HP = False
        enable_TESS = False
        enable_BESS = False

    elif energy_label == 'E':
        glazing_type = 'single'
        LRoof = 0.1
        LCavity_wall = 0.05
        Cavity_Filling = True
        Cavity_type = 'air'
        Wall_Cover = False
        
        enable_PV = False
        enable_HP = False
        enable_TESS = False
        enable_BESS = False

    elif energy_label == 'F':
        glazing_type = 'single'
        LRoof = 0.05
        LCavity_wall = 0.03
        Cavity_Filling = True
        Cavity_type = 'air'
        Wall_Cover = False 
        
        enable_PV = False
        enable_HP = False
        enable_TESS = False
        enable_BESS = False

    elif energy_label == 'G':
        glazing_type = 'single'
        LRoof = 0.03
        LCavity_wall = 0
        Cavity_Filling = False
        Cavity_type = 'air'
        Wall_Cover = False   
        
        enable_PV = False
        enable_HP = False
        enable_TESS = False
        enable_BESS = False
        
    return [glazing_type, LRoof, LCavity_wall, Cavity_Filling, Cavity_type, Wall_Cover, enable_PV, enable_HP, enable_TESS, enable_BESS]


def plot_sensitivity_analysis(registry, category_order, X_category, labels = ['House_type', 'Category', 'Thermal_demand', 'Gas_consumption', 'CAPEX', 'OPEX']):
    # from pandas import DataFrame
    import matplotlib.pyplot as plt
    
    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman"
    })
    
    House_type = [i[0] for i in registry]
    Category =  [i[-6] for i in registry]
    Thermal_demand = [i[-5] for i in registry]
    Gas_consumption = [i[-4] for i in registry]
    CAPEX = [i[-3] for i in registry]
    OPEX = [i[-2] for i in registry]
    
    # df = DataFrame([House_type, Category, Thermal_demand, Gas_consumption, CAPEX, OPEX]).T.set_axis(labels, axis=1)

    # plt.figure(constrained_layout=True)        
    # pivot_df = df.pivot(index='Category', columns='House_type', values='Thermal_demand')
    # pivot_df[category_order]
    # pivot_df.plot.bar()
    # plt.show()
    
    # return df
    

    from numpy import arange
    import matplotlib.pyplot as plt 
    
    barWidth = 0.25
    fig, ax = plt.subplots() 
    
    br1 = arange(len(category_order)) 
    br2 = [x + barWidth for x in br1] 
    br3 = [x + barWidth for x in br2] 
    
    bars1 = plt.bar(br1, [i for i,j in zip(Thermal_demand, House_type) if j == 'studio'], color ='r', width = barWidth, 
            edgecolor ='grey', label = 'Studio')
    bars2 = plt.bar(br2, [i for i,j in zip(Thermal_demand, House_type) if j == 'apartment'], color ='g', width = barWidth, 
            edgecolor ='grey', label = 'Apartment') 
    bars3 = plt.bar(br3, [i for i,j in zip(Thermal_demand, House_type) if j == 'stand_alone'], color ='b', width = barWidth, 
            edgecolor ='grey', label = 'Stand Alone') 

    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height, f'{height:.1f}', ha='center', va='bottom')

    
    plt.xlabel(X_category, fontsize = 15) 
    plt.ylabel('Thermal demand, $Q_{D}$, [MWh]', fontsize = 15) 
    plt.xticks([r + barWidth for r in range(len(category_order))], 
            category_order)
    
    plt.legend(ncol = 3, loc = 'lower center')
    plt.show()
    

    fig, ax = plt.subplots() 
    
    bars1 = plt.bar(br1, [i for i,j in zip(Gas_consumption, House_type) if j == 'studio'], color ='r', width = barWidth, 
            edgecolor ='grey', label = 'Studio')
    bars2 = plt.bar(br2, [i for i,j in zip(Gas_consumption, House_type) if j == 'apartment'], color ='g', width = barWidth, 
            edgecolor ='grey', label = 'Apartment') 
    bars3 = plt.bar(br3, [i for i,j in zip(Gas_consumption, House_type) if j == 'stand_alone'], color ='b', width = barWidth, 
            edgecolor ='grey', label = 'Stand Alone') 

    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height, f'{height:.1f}', ha='center', va='bottom')
    
    plt.xlabel(X_category, fontsize = 15) 
    plt.ylabel('Gas consumption, $V_{Gas}$, [m$^{3}$]', fontsize = 15) 
    plt.xticks([r + barWidth for r in range(len(category_order))], 
            category_order)
    
    plt.legend(ncol = 3, loc = 'lower center')
    plt.show()    
    

    fig, ax = plt.subplots() 
    
    bars1 = plt.bar(br1, [i for i,j in zip(CAPEX, House_type) if j == 'studio'], color ='r', width = barWidth, 
            edgecolor ='grey', label = 'Studio')
    bars2 = plt.bar(br2, [i for i,j in zip(CAPEX, House_type) if j == 'apartment'], color ='g', width = barWidth, 
            edgecolor ='grey', label = 'Apartment') 
    bars3 = plt.bar(br3, [i for i,j in zip(CAPEX, House_type) if j == 'stand_alone'], color ='b', width = barWidth, 
            edgecolor ='grey', label = 'Stand Alone') 

    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height, f'{height:.1f}', ha='center', va='bottom')
    
    plt.xlabel(X_category, fontsize = 15) 
    plt.ylabel('CAPEX, [€]', fontsize = 15) 
    plt.xticks([r + barWidth for r in range(len(category_order))], 
            category_order)
    
    plt.legend(ncol = 1, loc = 'upper right')
    plt.show()  


    fig, ax = plt.subplots() 
    
    bars1 = plt.bar(br1, [i for i,j in zip(OPEX, House_type) if j == 'studio'], color ='r', width = barWidth, 
            edgecolor ='grey', label = 'Studio')
    bars2 = plt.bar(br2, [i for i,j in zip(OPEX, House_type) if j == 'apartment'], color ='g', width = barWidth, 
            edgecolor ='grey', label = 'Apartment') 
    bars3 = plt.bar(br3, [i for i,j in zip(OPEX, House_type) if j == 'stand_alone'], color ='b', width = barWidth, 
            edgecolor ='grey', label = 'Stand Alone') 

    for bars in [bars1, bars2, bars3]:
        for bar in bars:
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height, f'{height:.1f}', ha='center', va='bottom')
    
    plt.xlabel(X_category, fontsize = 15) 
    plt.ylabel('OPEX, [€]', fontsize = 15) 
    plt.xticks([r + barWidth for r in range(len(category_order))], 
            category_order)
    
    plt.legend(ncol = 3, loc = 'lower right')
    plt.show() 
    
def preload_data(House_type, P_PV_av_studio, P_PV_av_apartment, P_PV_av_stand_alone, P_L_studio, P_L_apartment, P_L_stand_alone):
    
    if House_type == 'studio':
        return [P_PV_av_studio, P_L_studio]
        
    elif House_type == 'apartment':
        return [P_PV_av_apartment, P_L_apartment]
    
    elif House_type == 'stand_alone':
        return [P_PV_av_stand_alone, P_L_stand_alone]
    
    
#%%##############################################################################    

House_type = ['studio', 'apartment', 'stand_alone']    # 'studio', 'apartment', 'stand_alone'
enable_PV = False
enable_HP = False
enable_TESS = False
enable_BESS = False

glazing_type = ['single', 'double', 'triple']                   # 'single', 'double', 'triple'
LRoof = [0.03, 0.05, 0.1, 0.15, 0.2, 0.3]                 # 0.03, 0.05, 0.1, 0.15, 0.15, 0.2, 0.3

LCavity_wall = [0, 0.03, 0.05, 0.08, 0.1]           # 0, 0.03, 0.05, 0.05, 0.08, 0.08, 0.1
Cavity_Filling = [False, True]
Cavity_type = ['air', 'EPS']                                    # 'air', 'EPS'
Wall_Cover = [False, True]


Energy_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G']

windows_registry = []
roof_registry = []
wall_cavity_registry = []
wall_cover_registry = []
energy_labels_registry = []

reference_registry = []

# [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in] = simulate_house(5, House_type, enable_PV, enable_HP, enable_TESS, enable_BESS)


for house in House_type:
    
    print(house)
    # start_time = time.time()    # The timer is initializad.  
    # P_L = get_base_load(house)
    # [T_amb, T_set, T_soil, G, P_PV_av, df_Prices] = Load_data(P_L)
    [P_L, P_PV_av] = preload_data(house, P_PV_av_studio, P_PV_av_apartment, P_PV_av_stand_alone, P_L_studio, P_L_apartment, P_L_stand_alone)
    # end_time = time.time()    # The timer is finished
    # print('Loading time: ', (end_time - start_time)/60)
    
    
    # for window in glazing_type:
        
    #     [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in] = simulate_house(5, T_amb, T_set, T_soil, G, P_PV_av, P_L, df_Prices, house, enable_PV, enable_HP, enable_TESS, enable_BESS, glazing_type = window, Cavity_Filling = Cavity_Filling[0], Cavity_type = Cavity_type[0], Wall_Cover = Wall_Cover[0], LRoof = LRoof[3], LCavity_wall = LCavity_wall[2])
        
    #     windows_registry.append([house, window, Calculate_thermal_demand(Qdot_Boiler, Qdot_HP, Qdot_TESS), Calculate_gas_consumption(Qdot_Boiler), CAPEX(house, window, 0, False, 'air', False), OPEX(P_Grid, Qdot_Boiler, df_Prices), [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in]]) #
    
    # for roof in LRoof:
        
    #     [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in] = simulate_house(5, T_amb, T_set, T_soil, G, P_PV_av, P_L, df_Prices, house, enable_PV, enable_HP, enable_TESS, enable_BESS, glazing_type = glazing_type[1], Cavity_Filling = Cavity_Filling[0], Cavity_type = Cavity_type[0], Wall_Cover = Wall_Cover[0], LRoof = roof, LCavity_wall = LCavity_wall[2])
        
    #     roof_registry.append([house, roof, Calculate_thermal_demand(Qdot_Boiler, Qdot_HP, Qdot_TESS), Calculate_gas_consumption(Qdot_Boiler), CAPEX(house, 'single', roof, False, 'air', False), OPEX(P_Grid, Qdot_Boiler, df_Prices), [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in]]) # 
        
    # for cavity in Cavity_type:
    #     for lcavity in LCavity_wall:
    #         # print('Cavity_type main :', cavity, ' - LCavity_wall main: ', lcavity)
    #         [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in] = simulate_house(5, T_amb, T_set, T_soil, G, P_PV_av, P_L, df_Prices, house, enable_PV, enable_HP, enable_TESS, enable_BESS, glazing_type = glazing_type[1], Cavity_Filling = True, Cavity_type = cavity, Wall_Cover = Wall_Cover[0], LRoof = LRoof[3], LCavity_wall = lcavity)
        
    #         wall_cavity_registry.append([house, cavity, lcavity, Calculate_thermal_demand(Qdot_Boiler, Qdot_HP, Qdot_TESS), Calculate_gas_consumption(Qdot_Boiler), CAPEX(house, 'single', 0, True, cavity, False), OPEX(P_Grid, Qdot_Boiler, df_Prices), [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in]])
        
    # for cover in Wall_Cover:
        
    #     [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in] = simulate_house(5, T_amb, T_set, T_soil, G, P_PV_av, P_L, df_Prices, house, enable_PV, enable_HP, enable_TESS, enable_BESS, glazing_type = glazing_type[1], Cavity_Filling = Cavity_Filling[0], Cavity_type = Cavity_type[0], Wall_Cover = cover, LRoof = LRoof[3], LCavity_wall = LCavity_wall[2])
        
    #     wall_cover_registry.append([house, cover, Calculate_thermal_demand(Qdot_Boiler, Qdot_HP, Qdot_TESS), Calculate_gas_consumption(Qdot_Boiler), CAPEX(house, 'single', 0, False, 'air', cover), OPEX(P_Grid, Qdot_Boiler, df_Prices), [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in]]) #
        
    # for label in Energy_labels:        
    #     [window, roof, lcavity, filling, cavity, cover, enable_PV, enable_HP, enable_TESS, enable_BESS] = Energy_label_parameters(label)
    #     # [glazing_type, LRoof, LCavity_wall, Cavity_Filling, Cavity_type, Wall_Cover]
    #     # print('Energy label ', label, ' parameters  ', [window, roof, lcavity, filling, cavity, cover, enable_PV, enable_HP, enable_TESS, enable_BESS])
        
    #     [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in] = simulate_house(5, T_amb, T_set, T_soil, G, P_PV_av, P_L, df_Prices, house, enable_PV, enable_HP, enable_TESS, enable_BESS, glazing_type = window, Cavity_Filling = filling, Cavity_type = cavity, Wall_Cover = cover, LRoof = roof, LCavity_wall = lcavity)
        
    #     energy_labels_registry.append([house, label, Calculate_thermal_demand(Qdot_Boiler, Qdot_HP, Qdot_TESS), Calculate_gas_consumption(Qdot_Boiler), CAPEX(house, window, roof, filling, cavity, cover, enable_HP, max(P_PV_av)), OPEX(P_Grid, Qdot_Boiler, df_Prices), [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in]]) #


for house in House_type:
    
    print(house)
    # start_time = time.time()    # The timer is initializad.  
    # P_L = get_base_load(house)
    # [T_amb, T_set, T_soil, G, P_PV_av, df_Prices] = Load_data(P_L)
    [P_L, P_PV_av] = preload_data(house, P_PV_av_studio, P_PV_av_apartment, P_PV_av_stand_alone, P_L_studio, P_L_apartment, P_L_stand_alone)
    # end_time = time.time()    # The timer is finished
    # print('Loading time: ', (end_time - start_time)/60)

    
    [window, roof, lcavity, filling, cavity, cover, enable_PV, enable_HP, enable_TESS, enable_BESS] = Energy_label_parameters('D')
    enable_HP = True
    enable_PV = True
    [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in] = simulate_house(5, T_amb, T_set, T_soil, G, P_PV_av, P_L, df_Prices, house, enable_PV, enable_HP, enable_TESS, enable_BESS, glazing_type = window, Cavity_Filling = filling, Cavity_type = cavity, Wall_Cover = cover, LRoof = roof, LCavity_wall = lcavity)
    reference_registry.append([house, 'D+PV+HP', Calculate_thermal_demand(Qdot_Boiler, Qdot_HP, Qdot_TESS), Calculate_gas_consumption(Qdot_Boiler), CAPEX(house, window, roof, filling, cavity, cover, enable_HP, max(P_PV_av)), OPEX(P_Grid, Qdot_Boiler, df_Prices), [P_PV, P_BESS, P_HP, P_HP_TESS, P_Grid, SOC_BESS, C_BESS, Qdot_D, Qdot_HP, Qdot_HP_TESS, T_TESS, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, T_in]]) # 

# plot_sensitivity_analysis(windows_registry, glazing_type, 'Window type')
# plot_sensitivity_analysis(roof_registry, LRoof, 'Roof insolation')

# plot_sensitivity_analysis([i for i in wall_cavity_registry if i[1] == 'air'], LCavity_wall, 'Cavity wall insolation -  air')
# plot_sensitivity_analysis([i for i in wall_cavity_registry if i[1] == 'EPS'], LCavity_wall, 'Cavity wall insolation -  EPS')

# plot_sensitivity_analysis(wall_cover_registry, Wall_Cover, 'External wall insolation')

# plot_sensitivity_analysis(energy_labels_registry, Energy_labels, 'Energy Label')