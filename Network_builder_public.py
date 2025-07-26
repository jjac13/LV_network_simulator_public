# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a set of functions used to simulate a low-
## voltage residential network, that can or not be controlled, and which can
## include multi-carrier energy systems per house. The mathematical description
## of the library can be found in the paper:
## Alpízar-Castillo, J., Ramírez-Elizondo, L., van Voorden, A., and Bauer, P. 
## (2025) "Aggregated residential multi-carrier energy storage as voltage 
## control provider in low-voltage distribution networks", Journal of Energy
## Storage, 132.A, 117507.
## https://www.sciencedirect.com/science/article/pii/S2352152X25022200    
## Created by: Dr. Joel Alpízar Castillo.
## TU Delft
## Version: 2.1

###############################################################################
###########################   Imported modules   ##############################

import pandas as pd
import time
from numpy import arange, array

global t_registry
t_registry = []

###############################################################################
    
def Create_loads(DF_Network, Network_headers, Load_data = True, show_plots = False, csv_name = 'S_n_test.csv', csv_address = '', csv_delim = ',', t = 365*24*4):

    if Load_data == True:
        import csvreader
        S_n_Data = csvreader.read_data(csv = csv_name, address = csv_address, delim = csv_delim)
        S_n_Data.data2array()
        return S_n_Data.ar


    else:
        from re import findall
        from Probabilistic_Load_Profile_Generator import Create_profile
        from numpy import zeros
    
        S_n = zeros((DF_Network[Network_headers[1]].max(),t))
        for i, node, profile in zip(DF_Network[Network_headers[4]].isnull(), DF_Network[Network_headers[1]], DF_Network[Network_headers[4]]):
            if not i:
                S_n[node-1,:] = Create_profile(True, int(findall(r"\d+", profile)[0]))
    
        if show_plots:
            import matplotlib.pylab as plt
            plt.rcParams.update({
            #    "text.usetex": True,
                "font.family": "Times New Roman",
                'font.size': 16
            })
    
            
            plt.figure(constrained_layout=True)  
            plt.grid()
            for S, i, node, profile in zip(S_n, DF_Network[Network_headers[4]].isnull(), DF_Network[Network_headers[1]], DF_Network[Network_headers[4]]):
                if not i:
                    plt.plot(S, label = 'Node '+str(node)+', '+str(int(findall(r"\d+", profile)[0]))+' kWh')
            plt.xlim([0,t])
            plt.ylabel('Demanded power, $P_{Grid}$, [kW]')
            plt.xlabel('Time, t, [s]') 
            plt.legend()
            plt.show()      
        return S_n.transpose()            
    

###############################################################################   

def Create_admittance_matrix(DF_Network, Network_headers, show_heatmap = False):
    from numpy import zeros
    
    A = zeros((max(DF_Network[Network_headers[1]]), max(DF_Network[Network_headers[1]]))) 
    
    for i,j in zip(DF_Network[Network_headers[0]],DF_Network[Network_headers[1]]):

        A[j-1,i-1] = 1
        

    
    if show_heatmap: 

        import matplotlib.pylab as plt
        plt.rcParams.update({
            "font.family": "Times New Roman",
            'font.size': 16
        })

        
        plt.figure(constrained_layout=True)

        plt.imshow(A, cmap = 'binary')

        plt.title('Admittance Matrix')
        plt.xlabel('Node')
        plt.ylabel('Node')
        plt.show()        
    
    return A

###############################################################################

def Calculate_impedance(DF_Wires, distance, wire_label, Wires_headers, km = True):   
    
    if km:
        return distance*(DF_Wires[DF_Wires[Wires_headers[0]] == wire_label][Wires_headers[1]].item() + 1j*DF_Wires[DF_Wires[Wires_headers[0]] == wire_label][Wires_headers[2]].item())/1000
    
    else:
        return distance*(DF_Wires[DF_Wires[Wires_headers[0]] == wire_label][Wires_headers[1]].item() + 1j*DF_Wires[DF_Wires[Wires_headers[0]] == wire_label][Wires_headers[2]].item())

###############################################################################

def Create_impedance_matrix(DF_Network, Network_headers, Wires_headers, show_heatmap = False, wire_file = 'Wires.xlsx', km = True):
    from numpy import zeros
    
    DF_Wires = pd.read_excel(wire_file)
    
    
    Z = zeros((max(DF_Network[Network_headers[1]]), max(DF_Network[Network_headers[1]])),dtype=complex) 
     
    
    for i,j,k,l in zip(DF_Network[Network_headers[0]], DF_Network[Network_headers[1]], DF_Network[Network_headers[2]], DF_Network[Network_headers[3]]):        

        Z[j-1,i-1] = Calculate_impedance(DF_Wires, k, l, Wires_headers, km)


    if show_heatmap:
       
        import matplotlib.pylab as plt
        plt.rcParams.update({
        #    "text.usetex": True,
            "font.family": "Times New Roman",
            'font.size': 16
        }) 

        
        fig, (ax) = plt.subplots(1, 1, sharey=True) # , constrained_layout=True
        a = plt.imshow(abs(Z), cmap = 'binary')
        cbar = fig.colorbar(a, ax=ax)
        cbar.set_label('Impedance, [$\Omega$]')        
        
        plt.title('Impedance Matrix')
        plt.xlabel('Node')
        plt.ylabel('Node')
        plt.show()           
    
    return Z

###############################################################################

def plot_current(I):
    from numpy import arange
    import matplotlib.pylab as plt

    plt.rcParams.update({
        "font.family": "Times New Roman",
        'font.size': 16
    })  
    
    plt.figure(constrained_layout=True)
    plt.plot([abs(i)*i.real/abs(i.real) for i in I])
    plt.title('Current vector')
    plt.xlabel('Node')
    plt.xlim([1, len(I)-1])
    plt.xticks(arange(0, len(I), step=1), arange(1, len(I)+1, step=1))            
    plt.ylabel('Current [A]')
    plt.grid()
    plt.show()  

###############################################################################

def Estimate_Currents(A, S_n, I_0, V_0, dI = 0.01, show_plot = False):
    from numpy import divide, matmul, transpose, where, isnan
    

    
    if [] != where(isnan(divide(S_n, V_0))):      
        I_n = divide(S_n, V_0)
    else:
        print('zeros')
        I_n = array([S/V if V!=0 else 0 for S,V in zip(S_n, V_0)])


    I = matmul(transpose(A),transpose(I_0)) + I_n
    
    i = 0
    while max(abs(I - I_0)) > dI:
        
        I_0 = I
        I = matmul(transpose(A),I_0) + I_n
        
        if i>100:
            print('Error in current estimation')
            break

    if show_plot:
        plot_current(I)
        
    return I

###############################################################################

def plot_voltage(V_n, V_feeder = 400): 
    import matplotlib.pylab as plt

    plt.rcParams.update({
        "font.family": "Times New Roman",
        'font.size': 16
    })  
    
    plt.figure(constrained_layout=True)
    plt.plot([abs(i)/V_feeder for i in V_n])
    plt.title('Voltage vector')
    plt.xlabel('Node')
    plt.xlim([1, len(V_n)-1])    
    plt.ylabel('Voltage [pu]')
    plt.grid()
    plt.show()   

###############################################################################

def Estimate_Node_voltage(A, Z, S_n, I_0, V_0, V_feeder = 400, dV = 0.001, dI = 0.01, show_plot = False): # V_feeder = 400
    from numpy import matmul, zeros
    
    I_0 = Estimate_Currents(A, S_n, I_0, V_0, dI)
    
    B = zeros(len(I_0))
    B[0] = 1

    V_n = matmul(A,V_0) - matmul(Z,I_0) + B*V_feeder
    
    i = 0
    while max(abs(V_n - V_0)) > dV:
        
        V_0 = V_n
        I_0 = Estimate_Currents(A, S_n, I_0, V_0, dI)
        V_n = matmul(A,V_0) - matmul(Z,I_0) + B*V_feeder
        
        i+=1
        if i>1000:
            print('Error in voltage estimation: ', max(abs(V_n - V_0)))
            break    
    
    
    if show_plot:
      plot_voltage(V_n, V_feeder)
     
        
    return [V_n, I_0]

###############################################################################

def Create_network(DF_Network, Network_headers, Wires_headers, wire_file = 'Gaia_cables.xlsx'):
    A = Create_admittance_matrix(DF_Network, Network_headers, False)
    Z = Create_impedance_matrix(DF_Network, Network_headers, Wires_headers, False, wire_file)

    return [A, Z]    
    
###############################################################################

def Network_state(A, Z, df_S_n, timestep, I_0, V_0):
    
    S_n_timestep = df_S_n.loc[timestep].tolist()
    return Estimate_Node_voltage(A, Z, [1000*S for S in S_n_timestep], I_0, V_0)
    

###############################################################################
    
def Create_PV_DF(S_n = False, t_0 = 0, t_simulation = 24*4, H = 0, P_PV_ref = array([]), module_power_ref = 0.315, return_modules = False, dt = 0.25):
    import csvreader
    from pandas import DataFrame
    
    CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
    CSVDataPV.data2array()
    labels = [str(i+1) for i in range(len(S_n[0]))]
    

    # if P_PV_ref != False:#P_PV_ref.any():
    if P_PV_ref.any() != False:#P_PV_ref.any():        
        P_PV = DataFrame([[ref*i[0]/module_power_ref/1000 for i in CSVDataPV.ar[t_0:t_0+t_simulation+H]] for ref in P_PV_ref])    

    # elif P_PV_ref == False:#S_n.any():
    else:
        E_PV_ref_module = sum([i[0] for i in CSVDataPV.ar])*dt/1000

        
        n_modules = [abs(sum(S_n[:,node])*dt//E_PV_ref_module) for node in range(len(S_n[0]))]
        

        P_PV = DataFrame([[n*i[0]/1000 for i in CSVDataPV.ar[t_0:t_0+t_simulation+H]] for n in n_modules])
        
    
    P_PV = P_PV.T
    P_PV.columns = labels
    
    if not return_modules:
    
        return P_PV
    
    else: 
        return n_modules

###############################################################################

def BESS_degradation_LFP(SoC_0, SoC_f, Capacity_BESS, Capacity_BESS_BOL = 10, T = 20 + 273, V = 62.7, model = 'Wang', dt = 0.25):
    #Modelling the cycling degradation of Li-ion batteries: Chemistry influenced stress factors - 10.1016/j.est.2021.102765    
    from numpy import exp

    if model == 'Olmos':    
        SoC_0*=100
        SoC_f*=100
    
        DOD = abs(SoC_f - SoC_0)
        mSOC = (SoC_f + SoC_0)/2
        N = Capacity_BESS*DOD/(2*Capacity_BESS_BOL)
        
        if SoC_0 > SoC_f:   # Discharge
            Cch = 0
            Cdc = Capacity_BESS*(DOD/100)/(Capacity_BESS_BOL*dt)
    
        elif SoC_0 < SoC_f:   # Charge
            Cch = Capacity_BESS*(DOD/100)/(Capacity_BESS_BOL*dt)
            Cdc = 0
    
        else:               # Rest
            Cch = 0
            Cdc = 0       
            
        
        beta = 0.003414
        k_t = 5.8755
        k_dod = -0.0046
        k_ch = 0.1038
        k_dc = 0.296
        k_soc = 0.0513
        alpha = 0.869
        T_ref = 293
        soc_ref = 42
        
        delta = beta*exp(k_t*(T - T_ref)/T + k_dod*DOD + k_ch*Cch + k_dc*Cdc)*(1+k_soc*mSOC*(1-(mSOC/(2*soc_ref))))
        
        soh = 100 - delta*pow(N,alpha)
        return soh*Capacity_BESS/100


#########
    elif model == 'Vermeer':
        from math import sqrt
        
        isa = 1000*Capacity_BESS*abs(SoC_f - SoC_0)/(V*dt)
        
    
        c=[0.0008/3600, 0.39, 1.035, 50, 14.876/sqrt(24*3600)] # aging parameters
        R = 8.314 # Regnault constant 

        dt *= 3600
        
        ilossCycle = c[0]*c[2]/c[3]*exp(c[1]*abs(isa))*(1-SoC_0)*abs(isa)*dt # cyclic aging   
        ilossCal = c[4]*sqrt(dt)*exp(-24e3/R/T) # calendar aging
        iloss = ilossCycle + ilossCal # total aging
        
        dt /= 3600        
        
        return Capacity_BESS - iloss*V*dt
    
    elif model == 'Wang':
        from math import sqrt
        
        a = 8.61*1e-6
        b = -5.13*1e-3
        c = 7.63*1e-1
        d = -6.7*1e-3
        e = 2.35
        f = 14.876/sqrt(24)
        E_a = -24.5*1000
        R = 8.314 # Regnault constant         
        
        DOD = abs(SoC_f - SoC_0)        
        i_rate = DOD/dt # Capacity_BESS*(DOD)/(Capacity_BESS_BOL*dt)
        Ah_throughput = 1000*Capacity_BESS*DOD/(V) #1000*Capacity_BESS*DOD/(V*dt)
    
        
        ilossCycle = (a*T**2 + b*T + c) * exp((d*T+e)*i_rate) * Ah_throughput # cyclic aging   
        ilossCal = f * sqrt(dt) * exp(E_a/(R*T)) # calendar aging
        iloss = ilossCycle - ilossCal # total aging
        
        return Capacity_BESS*(100+iloss)/100

###############################################################################

def update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL = 10, charge_efficiency = 0.943, discharge_efficiency = 0.943, dt = 0.25):
    
    if P_BESS >= 0: # Discharge
        E_BESS = SoC_0*Capacity_BESS - dt*P_BESS/discharge_efficiency
    else:
        E_BESS = SoC_0*Capacity_BESS - dt*P_BESS*charge_efficiency
        
    
    return [E_BESS/Capacity_BESS, BESS_degradation_LFP(SoC_0, E_BESS/Capacity_BESS, Capacity_BESS, Capacity_BESS_BoL, model = 'Wang')]   


###############################################################################
################################   Heuristic  #################################

def Gas_Boiler(T_in, T_amb, T_set, mdot = 0.1, T_network = 50 + 273, efficiency = 0.8, c = 4200, House_type = 'apartment', glazing_type = 'double', Cavity_Filling = False, Cavity_type = 'air', Wall_Cover = False, LRoof = 0.2, LCavity_wall = 0.05):
    from MCES_library import new_house_Temperature, House_Thermal_Losses
    
    Qdot_D = House_Thermal_Losses(T_in, T_amb, House_type = House_type, glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[0]
    
    if T_in < T_set:
        Qdot_Boiler = mdot*c*(T_network - T_in)
    else:
        Qdot_Boiler = 0
    
    T_in_new = new_house_Temperature(T_in, House_Thermal_Losses(T_in, T_amb, House_type = House_type, glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[0], House_Thermal_Losses(T_in, T_amb, House_type = House_type, glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[1], Qdot_Boiler = Qdot_Boiler)

    return [Qdot_Boiler, T_in_new, Qdot_D]
###############################################################################
    


def Electric_Power_Balance(P_Load, P_HP, P_BESS, P_BESS_max, P_BESS_min, P_PV, P_grid_DSO, top = True, exact_match = True):

    
    if P_grid_DSO:
        if exact_match:
            if P_HP + P_Load - P_PV - P_BESS_max >= P_grid_DSO:         # The BESS and the PV cannot compensate (increase injection)
                return [P_PV, P_BESS_max, P_HP + P_Load - P_PV - P_BESS_max]

            elif (P_HP + P_Load - P_PV - P_BESS_max <= P_grid_DSO) and (P_HP + P_Load - P_PV >= P_grid_DSO):       # Curtailing the PV is not enough to compensate, but the BESS can compensate (increase injection)
                return [P_PV, P_HP + P_Load - P_PV - P_grid_DSO, P_grid_DSO]

            elif (P_HP + P_Load - P_PV <= P_grid_DSO) and (P_HP + P_Load >= P_grid_DSO):                    # Curtailing the PV is enough to compensate
                return [P_HP + P_Load - P_grid_DSO, 0, P_grid_DSO]

            elif P_HP + P_Load - P_BESS_min >= P_grid_DSO:              # The BESS has enough space to consume power to compensate (increase consumption)
                return [0, P_HP + P_Load - P_grid_DSO, P_grid_DSO]

            elif P_HP + P_Load - P_BESS_min <= P_grid_DSO:               # The BESS cannot consume enough power to compensate (increase consumption)          
                return [0, P_BESS_min, P_HP + P_Load - P_BESS_min]        
            
           
        else:
            if top:     # Peak-Shaving         
                if P_HP + P_Load - P_PV <= P_grid_DSO:      # No action required
                    return [P_PV, 0,  P_Load + P_HP - 0 - P_PV] 
                
                elif P_HP + P_Load - P_PV - P_BESS_max <= P_grid_DSO:   # BESS can compensate it
                    return [P_PV, P_HP + P_Load - P_PV - P_grid_DSO, P_grid_DSO]
                
                else:    # BESS cannot compensate completely
                    return [P_PV, P_BESS_max, P_Load + P_HP - P_BESS_max - P_PV]
        
            else:       # Power curtailment     
                if P_HP + P_Load - P_PV >= P_grid_DSO:     # No action required
                    return [P_PV, 0, P_Load + P_HP - 0 - P_PV]
                
                elif P_HP + P_Load >= P_grid_DSO:   # PV curtailment is enough
                    return [P_HP + P_Load - P_grid_DSO, 0, P_grid_DSO]
                
                elif P_HP + P_Load + P_BESS_min >= P_grid_DSO:    # BESS can compensate, curtailment is not enough
                    return [0, P_grid_DSO - P_HP - P_Load, P_grid_DSO]            
                
                else:   # BESS cannot compensate
                    return [0, P_BESS_min, P_Load + P_HP - P_BESS_min - P_PV]            
                
    else:           # No request from the DSO   
        return [P_PV, P_BESS, P_Load + P_HP - 0 - P_PV]     


###############################################################################

def check_P_BESS_setpoint(SoC_0, Capacity_BESS, P_Load, P_HP, P_HP_TESS, P_PV_av, P_grid_DSO, P_BESS_max = 10, dt = 0.25):
    from MCES_library import BESS_perm_max, BESS_perm_min    
    
    P_BESS_min = BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)
    P_BESS_max = BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)
    P_BESS_set_min = P_Load + P_HP + P_HP_TESS - P_PV_av - P_grid_DSO
    P_BESS_set_max = P_Load + P_HP + P_HP_TESS - 0 - P_grid_DSO
    
    if (P_BESS_set_min >= P_BESS_min) and (P_BESS_set_min <= P_BESS_max): # and (P_BESS_set_max <= P_BESS_max):
        return [P_BESS_set_min, P_PV_av, True]
    
    elif (P_BESS_set_min <= P_BESS_min) and (P_BESS_set_max >= P_BESS_min):
        return [P_BESS_min, P_BESS_min - P_BESS_set_min, True]
    
    elif (P_BESS_set_min >= P_BESS_max):
        return [P_BESS_max, 0, False]
    
    else: # (P_BESS_set_max <= P_BESS_min)
        return [P_BESS_min, P_PV_av, False]
    
###############################################################################    

def house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO = False, top = True, P_BESS_ref = 0, Capacity_BESS_BoL = 10, P_BESS_max = 10, T_TESS_max = 90 + 273, setpoint = False):
    from MCES_library import BESS_perm_max, BESS_perm_min, HP_Power, update_TESS

    if setpoint:
        [P_PV, P_BESS, P_Grid] = Electric_Power_Balance(P_Load, HP_Power(HP_active, T_amb, T_ret)[0] + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_BESS_ref, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), P_PV_av, P_grid_DSO, top)
    
    else:
        [P_PV, P_BESS, P_Grid] = [P_PV_av, P_BESS_ref, P_Load + HP_Power(HP_active, T_amb, T_ret)[0] + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0] - P_PV_av - P_BESS_ref]
    
    [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1])                               
    [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)

    return [HP_Power(HP_active, T_amb, T_ret)[0], HP_Power(HP_active, T_amb, T_ret)[1], P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS]
    
###############################################################################    
def keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = 16+273):
    if T_in > T_in_min:
        return [TESS_active, HP_active, HP_charge_TESS]
    else:
        return [False, True, False]

###############################################################################
                        
def heuristic_control(T_in, T_set, T_amb, T_TESS_0, T_ret, T_soil, G, P_Load, P_grid_DSO, top, SoC_0, Cost_grid, Cost_grid_median, Cost_grid_Q1, Cost_grid_Q3, Capacity_BESS, P_PV_av, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP = True, enable_TESS = True, enable_HP_2_TESS = True, external_control = False, SoC_BESS_min = 0.9, SoC_BESS_max = 0.2, m = 4000, TESS_min_op_T = 55 + 273, TESS_max_op_T = 90 + 273, T_TESS_min = 50 + 273, T_TESS_max = 90 + 273, T_in_crit = 16 + 273, m_dot = 0.22, c_f = 4200, dt = 0.25, House_type = 'apartment', force_temperature = False, glazing_type = 'double', Cavity_Filling = False, Cavity_type = False, Wall_Cover = False, LRoof = 0.2, LCavity_wall = 0.05):
    from MCES_library import BESS_perm_max, BESS_perm_min, HP_Power, new_house_Temperature, House_Thermal_Losses # , update_TESS
    
    if not force_temperature:
        T_in_min = 16+273 
    else:
        T_in_min = T_set
    
    if not external_control:
        # print('base')
        if T_in < T_set:    # Heating required
            if Cost_grid <= 0: # Negative prices, charge BESS
                TESS_active = enable_TESS*False
                HP_active = enable_HP*True
                HP_charge_TESS = enable_TESS*enable_HP*(not HP_active)*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                    
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, 0, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), Capacity_BESS_BoL, P_BESS_max, T_TESS_max)
               # print('No external control - Heating required - Negative prices, charge BESS')

    
            elif Cost_grid <= Cost_grid_Q1:     # Low prices, charge BESS
                TESS_active = enable_TESS*False
                HP_active = enable_HP*True
                HP_charge_TESS = enable_TESS*enable_HP*(not HP_active)*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4, Capacity_BESS_BoL, P_BESS_max, T_TESS_max)
               # print('No external control - Heating required - Low prices, charge BESS')
            
    
            elif Cost_grid <= Cost_grid_Q3:     # Average prices, BESS compensate
                TESS_active = enable_TESS*True
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*False*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_BESS_ref, P_PV, setpoint] = check_P_BESS_setpoint(SoC_0, Capacity_BESS, P_Load, HP_Power(HP_active, T_amb, T_ret)[0], HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_PV_av, P_grid_DSO, P_BESS_max = P_BESS_max)
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, P_BESS_ref, Capacity_BESS_BoL, P_BESS_max, T_TESS_max, setpoint = True)        

                
            else:                                   # High prices, discharge BESS
                TESS_active = enable_TESS*True
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*False*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4, Capacity_BESS_BoL, P_BESS_max, T_TESS_max)                
               # print('No external control - Heating required - High prices, discharge BESS')
                
            
        else:               # Heating not required
            if (Cost_grid <= 0)*(not external_control): # Negative prices, charge BESS
                TESS_active = enable_TESS*False
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*(T_TESS_0 < T_TESS_max)*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, 0, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), Capacity_BESS_BoL, P_BESS_max, T_TESS_max)                
                # print('No external control - Heating not required - Negative prices, charge BESS')

    
            elif (Cost_grid <= Cost_grid_Q1)*(not external_control):     # Low prices, charge BESS
                TESS_active = enable_TESS*False
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*(T_TESS_0 < T_TESS_max)*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4, Capacity_BESS_BoL, P_BESS_max, T_TESS_max)                
#                print('No external control - Heating not required - Low prices, charge BESS')
    
            elif Cost_grid <= Cost_grid_Q3:     # Average prices, BESS compensate
                TESS_active = enable_TESS*False
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*False*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_BESS_ref, P_PV, setpoint] = check_P_BESS_setpoint(SoC_0, Capacity_BESS, P_Load, HP_Power(HP_active, T_amb, T_ret)[0], HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_PV_av, P_grid_DSO, P_BESS_max = P_BESS_max)
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, P_BESS_ref, Capacity_BESS_BoL, P_BESS_max, T_TESS_max, setpoint = True)        
#                print('No external control - Heating not required - Average prices, BESS compensate')
                
            else:                                   # High prices, discharge BESS
                TESS_active = enable_TESS*False
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*False*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4, Capacity_BESS_BoL, P_BESS_max, T_TESS_max)
#                print('No external control - Heating not required - High prices, discharge BESS')

    else:   # External control
#################

        # General:
        if Cost_grid <= Cost_grid_Q1:
            delta_HP_e = 1
            delta_TESS_e = 0
            delta_HP_TESS_e = 1        
        else:
            delta_HP_e = 0
            delta_TESS_e = 1
            delta_HP_TESS_e = 0
        
        if T_in <= T_set:
            delta_HP_T = 1
            delta_TESS_T = 1
        else:
            delta_HP_T = 0
            delta_TESS_T = 0
            
        delta_HP = delta_HP_e*delta_HP_T
        delta_TESS = delta_TESS_e*delta_TESS_T

        if T_TESS_0 < T_TESS_max:
            delta_HP_TESS_t = 1
        else:
            delta_HP_TESS_t = 0
            
        if delta_HP == 1:
            delta_HP_TESS_HP = 0
        else:
            delta_HP_TESS_HP = 1
            
        delta_HP_TESS = delta_HP_TESS_e*delta_HP_TESS_t*delta_HP_TESS_HP
        
        TESS_active = enable_TESS*delta_TESS
        HP_active = enable_HP*delta_HP
        HP_charge_TESS = enable_TESS*enable_HP*enable_HP_2_TESS*delta_HP_TESS
        
        [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
        
        if external_control: # external control
            [P_BESS_ref, P_PV, setpoint] = check_P_BESS_setpoint(SoC_0, Capacity_BESS, P_Load, HP_Power(HP_active, T_amb, T_ret)[0], HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_PV_av, P_grid_DSO, P_BESS_max = P_BESS_max)
            
            if not setpoint and T_in > 16+273:
                TESS_active = False
                HP_active = False
                HP_charge_TESS = False
            
                [P_BESS_ref, P_PV, setpoint] = check_P_BESS_setpoint(SoC_0, Capacity_BESS, P_Load, HP_Power(HP_active, T_amb, T_ret)[0], HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_PV_av, P_grid_DSO, P_BESS_max = P_BESS_max)
        
        if not external_control: # not external control
            # print('Controlled')
            if Cost_grid <= 0:
                P_PV = 0
                P_BESS_ref = BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)
            elif Cost_grid <= Cost_grid_Q1:
                P_PV = P_PV_av
                P_BESS_ref = BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4
            elif Cost_grid <= Cost_grid_Q3:
                P_PV = P_PV_av
                P_BESS_ref = 0
            else:
                P_PV = P_PV_av
                P_BESS_ref = BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4
                    
                
        
        [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, P_BESS_ref, Capacity_BESS_BoL, P_BESS_max, T_TESS_max, setpoint = True)        
        
    if enable_HP:    
        Qdot_Boiler = 0
    else:
        Qdot_Boiler = Gas_Boiler(T_in, T_amb, T_set, House_type = House_type, glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[0]
    
    T_in = new_house_Temperature(T_in, House_Thermal_Losses(T_in, T_amb, House_type = House_type, glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[0], House_Thermal_Losses(T_in, T_amb, House_type = House_type, glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[1], Qdot_HP = Qdot_HP, Qdot_TESS = Qdot_TESS, Qdot_Boiler = Qdot_Boiler)
    
    if T_TESS_new < TESS_min_op_T:                
        enable_HP_2_TESS = True
    elif T_TESS_new > TESS_max_op_T:
        enable_HP_2_TESS = False  
    

    return [P_PV, P_BESS, P_HP, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_Grid, SoC_BESS, Capacity_BESS, House_Thermal_Losses(T_in, T_amb, House_type = House_type, glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[0], Qdot_HP, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1], T_TESS_new, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, enable_HP_2_TESS, T_in]



###############################################################################

def GA_control(T_in, T_set, T_amb, T_TESS_0, T_ret, T_soil, G, P_Load, P_grid_DSO, top, SoC_0, costs, Capacity_BESS, P_PV_av, horizon = 0, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP_2_TESS = True, external_control = False, SoC_BESS_min = 0.9, SoC_BESS_max = 0.2, m = 4000, TESS_min_op_T = 55 + 273, TESS_max_op_T = 90 + 273, T_TESS_min = 50 + 273, T_TESS_max = 90 + 273, m_dot = 0.22, c_f = 4200, dt = 0.25, House_type = 'apartment'):
    from MCES_library import HP_Power, update_TESS, new_house_Temperature, House_Thermal_Losses
    from GA_EMS import GA_Optimization
       
    best_candidate = GA_Optimization(T_amb, T_in, T_set, T_soil, P_PV_av, P_Load, G, T_TESS_0, SoC_0, horizon, costs) # consecutive_generations = 5, individuals = 200, beta = 6, theta_E = 0.483*0.25/2, theta_T = 1/1, theta_CO2 = 1/2.5)
    
    
    [HP_charge_TESS, TESS_active, HP_active, P_PV, P_BESS] = [best_candidate[0], best_candidate[1], best_candidate[2], best_candidate[3], best_candidate[4]]
    [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
    [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb[0], T_ret)
    [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb[0], T_TESS_0, T_sup = T_TESS_max)[1])
    T_in = new_house_Temperature(T_in, House_Thermal_Losses(T_in, T_amb[0], House_type = House_type )[0], House_Thermal_Losses(T_in, T_amb[0], House_type = House_type )[1], Qdot_HP = Qdot_HP, Qdot_TESS = Qdot_TESS)
    
    
    return [P_PV, P_BESS, P_HP, HP_Power(HP_charge_TESS, T_amb[0], T_TESS_0, T_sup = T_TESS_max)[0], P_Load[0] + P_HP + HP_Power(HP_charge_TESS, T_amb[0], T_TESS_0, T_sup = T_TESS_max)[0] + P_BESS - P_PV, SoC_BESS, Capacity_BESS, House_Thermal_Losses(T_in, T_amb[0], House_type = House_type)[0], Qdot_HP, HP_Power(HP_charge_TESS, T_amb[0], T_TESS_0, T_sup = T_TESS_max)[1], T_TESS_new, Qdot_TESS, Qdot_TESS_SD, enable_HP_2_TESS, T_in]
    


###############################################################################


def Network_Control(A, Z, I_0, V_0, Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control, P_BESS_max, Capacity_BESS_BoL, enable_HP, enable_TESS, V_ref_min = 0.95*400, V_ref_max = 1.05*400, P_ref_min = -1.48, P_ref_max = 1.38, force_temperature = False):
    from numpy import zeros
    B = zeros(len(I_0))
    B[0] = 1  
        
    [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in] = Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control = False, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, enable_HP = enable_HP, enable_TESS = enable_TESS, force_temperature = force_temperature)

    if not follow_control:
#        print('Not controlled')
        return [[df_P_ref, df_DSO_operation], [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in]]
    
    else: 
#        print('Controlled')
        [V_n, I] = Network_state(A, Z, df_P_Grid, timestep, I_0, V_0)
        
        if (min(V_n) > V_ref_min) and (max(V_n) > V_ref_min):   # No action required
            return [[df_P_ref, df_DSO_operation], [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in]]

        else:
            P_ref = Find_power_setpoints(A, Z, B, V_n, I, df_P_Grid.loc[timestep].tolist(), selection+remaining, selection)
            for node in selection:
                df_DSO_operation.iloc[timestep,node] = True
                df_P_ref.iloc[timestep,node] = P_ref[node-1]
            return [[df_P_ref, df_DSO_operation], Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, enable_HP = enable_HP, enable_TESS = enable_TESS, force_temperature = force_temperature)]        


###############################################################################
def Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, ts, control_type = 'heuristic', follow_control = False, horizon = 0, energy_costs = [0, 0, 0, 0, 0.13, 0.483], CO2_costs = [0, 0, 0, 0, 0, 0.325], T_ret = 38 + 273, method = False, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP = True, enable_TESS = True, force_temperature = False, glazing_type = 'double', Cavity_Filling = False, Cavity_type = 'air', Wall_Cover = False, LRoof = 0.2, LCavity_wall = 0.05):
    
    if control_type == 'heuristic':
        if follow_control:
#            print('--- follow control ----')
            for node in range(len(df_T_in.loc[0])):
                if node in selection:
           
                    [df_P_PV.iloc[ts, node], df_P_BESS.iloc[ts, node], df_P_HP.iloc[ts, node], df_P_HP_TESS.iloc[ts, node], df_P_Grid.iloc[ts, node], df_SOC_BESS.iloc[ts, node], df_C_BESS.iloc[ts, node], df_Qdot_D.iloc[ts, node], df_Qdot_HP.iloc[ts, node], df_Qdot_HP_TESS.iloc[ts, node], df_T_TESS.iloc[ts, node], df_Qdot_TESS.iloc[ts, node], df_Qdot_TESS_SD.iloc[ts, node], df_Qdot_Boiler.iloc[ts, node], df_enable_HP_2_TESS.iloc[ts, node], df_T_in.iloc[ts, node]] = heuristic_control(df_T_in.iloc[ts-1, node], T_set[ts], T_amb[ts], df_T_TESS.iloc[ts-1, node], 38+273, T_soil[ts], 0, df_P_L.iloc[ts, node], df_P_ref.iloc[ts, node], df_DSO_operation.iloc[ts, node], df_SOC_BESS.iloc[ts-1, node], df_Prices['Current'][ts], df_Prices['Median'][ts], df_Prices['Q25'][ts], df_Prices['Q75'][ts], df_C_BESS.iloc[ts-1, node], df_P_PV_av.iloc[ts, node], enable_HP = enable_HP, enable_TESS = enable_TESS, enable_HP_2_TESS = df_enable_HP_2_TESS.iloc[ts-1, node], external_control = follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, House_type = df_House_type[Network_headers[5]][node+1], force_temperature = force_temperature)
                elif node in remaining:
                    [df_P_PV.iloc[ts, node], df_P_BESS.iloc[ts, node], df_P_HP.iloc[ts, node], df_P_HP_TESS.iloc[ts, node], df_P_Grid.iloc[ts, node], df_SOC_BESS.iloc[ts, node], df_C_BESS.iloc[ts, node], df_Qdot_D.iloc[ts, node], df_Qdot_HP.iloc[ts, node], df_Qdot_HP_TESS.iloc[ts, node], df_T_TESS.iloc[ts, node], df_Qdot_TESS.iloc[ts, node], df_Qdot_TESS_SD.iloc[ts, node], df_Qdot_Boiler.iloc[ts, node], df_enable_HP_2_TESS.iloc[ts, node], df_T_in.iloc[ts, node]] = [0, 0, 0, 0, df_P_L.iloc[ts, node], 0, 0, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1])[2], 0, 0, 0, 0, 0, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[0], False, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[1]]
                    
        else:
            for node in range(len(df_T_in.loc[0])):
                if node in selection:
           
                    [df_P_PV.iloc[ts, node], df_P_BESS.iloc[ts, node], df_P_HP.iloc[ts, node], df_P_HP_TESS.iloc[ts, node], df_P_Grid.iloc[ts, node], df_SOC_BESS.iloc[ts, node], df_C_BESS.iloc[ts, node], df_Qdot_D.iloc[ts, node], df_Qdot_HP.iloc[ts, node], df_Qdot_HP_TESS.iloc[ts, node], df_T_TESS.iloc[ts, node], df_Qdot_TESS.iloc[ts, node], df_Qdot_TESS_SD.iloc[ts, node], df_Qdot_Boiler.iloc[ts, node], df_enable_HP_2_TESS.iloc[ts, node], df_T_in.iloc[ts, node]] = heuristic_control(df_T_in.iloc[ts-1, node], T_set[ts], T_amb[ts], df_T_TESS.iloc[ts-1, node], 38+273, T_soil[ts], 0, df_P_L.iloc[ts, node], df_P_ref.iloc[ts, node], df_DSO_operation.iloc[ts, node], df_SOC_BESS.iloc[ts-1, node], df_Prices['Current'][ts], df_Prices['Median'][ts], df_Prices['Q25'][ts], df_Prices['Q75'][ts], df_C_BESS.iloc[ts-1, node], df_P_PV_av.iloc[ts, node], enable_HP = enable_HP, enable_TESS = enable_TESS, enable_HP_2_TESS = df_enable_HP_2_TESS.iloc[ts-1, node], external_control = follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, House_type = df_House_type[Network_headers[5]][node+1], force_temperature = force_temperature)
    
                elif node in remaining:
                    [df_P_PV.iloc[ts, node], df_P_BESS.iloc[ts, node], df_P_HP.iloc[ts, node], df_P_HP_TESS.iloc[ts, node], df_P_Grid.iloc[ts, node], df_SOC_BESS.iloc[ts, node], df_C_BESS.iloc[ts, node], df_Qdot_D.iloc[ts, node], df_Qdot_HP.iloc[ts, node], df_Qdot_HP_TESS.iloc[ts, node], df_T_TESS.iloc[ts, node], df_Qdot_TESS.iloc[ts, node], df_Qdot_TESS_SD.iloc[ts, node], df_Qdot_Boiler.iloc[ts, node], df_enable_HP_2_TESS.iloc[ts, node], df_T_in.iloc[ts, node]] = [0, 0, 0, 0, df_P_L.iloc[ts, node], 0, 0, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1])[2], 0, 0, 0, 0, 0, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[0], False, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[1]]


    elif control_type == 'GA':
        from GA_EMS import TS_forecast
        
        h = ts + horizon + 1        
        
        horizon_costs = [[],[]]
        for i in range(horizon + 1):
            energy_costs[-1] = df_Prices['Current'][ts + horizon]
            horizon_costs[0].append(energy_costs)        
            horizon_costs[1].append(CO2_costs)
            
            
            for node in range(len(df_T_in.loc[0])):            
                if node in selection:    
                    

                    [df_P_PV.iloc[ts, node], df_P_BESS.iloc[ts, node], df_P_HP.iloc[ts, node], df_P_HP_TESS.iloc[ts, node], df_P_Grid.iloc[ts, node], df_SOC_BESS.iloc[ts, node], df_C_BESS.iloc[ts, node], df_Qdot_D.iloc[ts, node], df_Qdot_HP.iloc[ts, node], df_Qdot_HP_TESS.iloc[ts, node], df_T_TESS.iloc[ts, node], df_Qdot_TESS.iloc[ts, node], df_Qdot_TESS_SD.iloc[ts, node], df_enable_HP_2_TESS.iloc[ts, node], df_T_in.iloc[ts, node]] = GA_control(
                        T_in = df_T_in.iloc[ts-1, node], 
                        T_set = T_set[ts:h], 
                        T_amb = TS_forecast(T_amb[ts:h], ts, method, forecast_type = 'Temperature'), 
                        T_TESS_0 = df_T_TESS.iloc[ts-1, node],
                        T_ret = T_ret,
                        T_soil = [min(i)+273 for i in T_soil[ts:h]],
                        P_PV_av = TS_forecast(df_P_PV_av.iloc[ts:h, node].tolist(), ts, method, forecast_type = 'PV'),
                        G = TS_forecast(G[ts:h], ts, method, forecast_type = 'Radiation'),
                        P_Load = TS_forecast(df_P_L.iloc[ts:h, node].tolist(), ts, method, forecast_type = 'Load'),
                        P_grid_DSO = df_P_ref.iloc[ts, node],
                        top = df_DSO_operation.iloc[ts, node],
                        SoC_0 = df_SOC_BESS.iloc[ts-1, node],
                        costs = horizon_costs,
                        Capacity_BESS = df_C_BESS.iloc[ts-1, node])

                elif node in remaining:
                    [df_P_PV.iloc[ts, node], df_P_BESS.iloc[ts, node], df_P_HP.iloc[ts, node], df_P_HP_TESS.iloc[ts, node], df_P_Grid.iloc[ts, node], df_SOC_BESS.iloc[ts, node], df_C_BESS.iloc[ts, node], df_Qdot_D.iloc[ts, node], df_Qdot_HP.iloc[ts, node], df_Qdot_HP_TESS.iloc[ts, node], df_T_TESS.iloc[ts, node], df_Qdot_TESS.iloc[ts, node], df_Qdot_TESS_SD.iloc[ts, node], df_Qdot_Boiler.iloc[ts, node], df_enable_HP_2_TESS.iloc[ts, node], df_T_in.iloc[ts, node]] = [0, 0, 0, 0, df_P_L.iloc[ts, node], 0, 0, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[2], 0, 0, 0, 0, 0, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[0], False, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[1]]

    return [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in]

###############################################################################

def assign_DF_values(DF, row, columns, value = 0):
    
    for col in range(len(DF.iloc[row])):
        if col in columns:
            DF.iloc[row, col] = value
    
    return DF
    
###############################################################################
    
def set_initialConditions(Network_headers, selection, remaining, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_T_in, T_amb, T_soil, df_House_type, initial_conditions = [0, 0, 0, 0, True, 0.5, 10, 0, 0,  0, True, 75+273, 19+273]):
    from MCES_library import House_Thermal_Losses, Qdot_SD_TESS2Soil
    
    [df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_T_in] = [assign_DF_values(DF, 0, selection, value) for DF,value in zip([df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_T_in], initial_conditions)]
        
    df_P_Grid.loc[0] = df_P_L.loc[0]
    df_P_ref = assign_DF_values(df_P_ref, 0, selection, 1.5)
    
    for node in selection:
        df_Qdot_D.iloc[0, node] = House_Thermal_Losses(df_T_in.iloc[0, node], T_amb[0], House_type = df_House_type[Network_headers[5]][node+1])[0]
        df_Qdot_TESS_SD.iloc[0, node] = Qdot_SD_TESS2Soil(T_soil[0], initial_conditions[10] - 273)

    for node in remaining:
        df_T_in.iloc[0, node] = 19+273
        df_Qdot_D.iloc[0, node] = House_Thermal_Losses(df_T_in.iloc[0, node], T_amb[0], House_type = df_House_type[Network_headers[5]][node+1])[0]
    
    return [df_P_HP, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_T_in]

###############################################################################

def create_pricesDF(Energy_price, t_simulation, t_0 = 0, dt = 0.25):
    from pandas import DataFrame
    from statistics import median
    from numpy import empty, percentile
    from math import floor

    labels = ['Min', "Q25", "Median", "Current", "Q75", "Max"]
    df_Prices = DataFrame(empty([int(t_simulation), len(labels)]), columns = labels) 

    for step in range(t_0, t_0 + t_simulation):       
        df_Prices['Min'][step-t_0] = min(Energy_price[int(floor(step/(24/dt))*24/dt):int((floor(step/(24/dt))+1)*24/dt)])
        df_Prices['Q25'][step-t_0] = percentile(Energy_price[int(floor(step/(24/dt))*24/dt):int((floor(step/(24/dt))+1)*24/dt)], 25)
        df_Prices['Median'][step-t_0] = median(Energy_price[int(floor(step/(24/dt))*24/dt):int((floor(step/(24/dt))+1)*24/dt)])    
        df_Prices['Current'][step-t_0] = Energy_price[step]
        df_Prices['Q75'][step-t_0] = percentile(Energy_price[int(floor(step/(24/dt))*24/dt):int((floor(step/(24/dt))+1)*24/dt)], 75)    
        df_Prices['Max'][step-t_0] = max(Energy_price[int(floor(step/(24/dt))*24/dt):int((floor(step/(24/dt))+1)*24/dt)])
        
    return df_Prices
        
    

###############################################################################    
def simulate_Network(DF_Network, Network_headers, Wires_headers, selection, remaining, S_n, t_simulation, V_0 = False, I_0 = False, t_0 = 0, wire_file = 'Wires.xlsx', H = 0, control_type = 'heuristic', follow_control = False, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP = True, enable_TESS = True, V_feeder = 400, dt = 0.25, force_temperature = False, start_day = 0, end_day = 1):
    from numpy import arange, array, empty, zeros, full
    from pandas import DataFrame
    import csvreader
    
#    print('Creating network')    
    
    [A, Z] = Create_network(DF_Network, Network_headers, Wires_headers, wire_file)

    
#    print('Loading data')
#    CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
    CSVDataTamb = csvreader.read_data(csv='Tamb_15min.csv', address='')
#    CSVDataP_Load = csvreader.read_data(csv='Load_Profile_15min.csv', address='', delim=',')
    CSVDataRad = csvreader.read_data(csv='Radiation_1min.csv', address='')
    CSVDataPrices = csvreader.read_data(csv='DA_Prices_15min.csv', address='')
    CSVDataTsoil = csvreader.read_data(csv='Soil_dy_10cm.csv', address='')
#    CSVDataPV.data2array()
    CSVDataTamb.data2array()
#    CSVDataP_Load.data2array()
    CSVDataRad.data2array()
    CSVDataPrices.data2array()
    CSVDataTsoil.data2cols()
#    P_PV_av = [i[0]*n_modules*module_power/module_power_ref/1000 for i in CSVDataPV.ar]
    T_amb = [i[0]+273 for i in CSVDataTamb.ar[t_0:t_0+t_simulation+H]]
#    P_Load = [i for i in CSVDataP_Load.ar[0][t_0:t_0+t_simulation+H]]
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

#    print('Creating dataframes')   
    
    labels = [str(i+1) for i in range(len(S_n[0]))]

    # Only for nodes with load
    df_P_PV_av = Create_PV_DF(S_n, t_0, t_simulation)
    df_P_PV = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)    
    df_P_L = DataFrame(S_n[t_0:t_0+t_simulation+H], columns = labels)#.div(1000)
    df_P_HP = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_P_HP_TESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_P_BESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_P_Grid = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_P_ref = DataFrame(full((int(t_simulation), len(A)), 1.4), columns = labels)
    df_DSO_operation = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)

    df_SOC_BESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_C_BESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)

    df_Qdot_D = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_Qdot_HP = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_Qdot_HP_TESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_Qdot_TESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_Qdot_TESS_SD = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_Qdot_Boiler = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    
    df_enable_HP_2_TESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)

    df_T_TESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    
    df_T_in = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    
    # For every node
    df_P_n = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_V_n = DataFrame(zeros([int(t_simulation), len(A)]))
    df_I_n = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    
    df_Prices = create_pricesDF(Energy_price, t_simulation)
    
    df_House_type = DF_Network[[Network_headers[1], Network_headers[5]]].sort_values(by=[Network_headers[1]])
    df_House_type = df_House_type.set_index(Network_headers[1])
    
    [df_P_HP, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_T_in] = set_initialConditions(Network_headers, selection, remaining, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_T_in, T_amb, T_soil, df_House_type)
    
    V_n_registry = empty([len(A), int(t_simulation)])
    I_registry = empty([len(A), int(t_simulation)]) 

   
    if not V_0:
        V_0 = zeros(len(A),dtype=complex) + V_feeder
        
    if not I_0:
        I_0 = zeros(len(A),dtype=complex)
        I_0 = Estimate_Currents(A, [1000*S for S in df_P_Grid.loc[0].tolist()], I_0, V_0)        
            
    [V_n, I] = Network_state(A, Z, df_P_Grid, 0, I_0, V_0) # + t_0
        
    V_n_registry[:,0] = V_n
    I_registry[:,0] = I    
    
    P_ref_min = -1.5
    P_ref_max = 1.5
        
    # print('Start Simulation')           

    for timestep in range(1,int(t_simulation)):
#        print('---', timestep, '---') # , 4*df_Prices['Current'][timestep]
        if (timestep)%(24/dt) == 0:
            print('Day', t_0/(24/dt) + int(timestep/(24/dt)))
        start_ts = time.time()
        
        [[df_P_ref, df_DSO_operation], [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in]] = Network_Control(A, Z, I, V_n, Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, enable_HP = enable_HP, enable_TESS = enable_TESS, P_ref_min = P_ref_min, P_ref_max = P_ref_max, force_temperature = force_temperature)
        
        [V_n, I] = Network_state(A, Z, df_P_Grid, timestep, I_0, V_0)
        
        V_n_registry[:,timestep] = V_n
        I_registry[:,timestep] = I
        
        
        end_ts = time.time()
        t_registry.append(end_ts-start_ts)

    df_V_n = DataFrame(V_n_registry).T
    df_V_n.columns = labels
    df_I_n = DataFrame(I_registry).T
    df_I_n.columns = labels

    return [V_n_registry, I_registry, [df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in, df_P_n, df_V_n, df_I_n, df_Prices, df_House_type]]

###############################################################################

def simulate_Network_Centralized_BESS(DF_Network, Network_headers, Wires_headers, selection, remaining, S_n, P_PV_peak, df_P_L, t_simulation, V_0 = False, I_0 = False, t_0 = 0, wire_file = 'Wires.xlsx', H = 0, control_type = 'heuristic', follow_control = False, P_BESS_max = 1000, Capacity_BESS_BoL = 1000, V_feeder = 400, dt = 0.25, V_ref_min = 0.95*400, V_ref_max = 1.05*400):
    from numpy import arange, array, empty, zeros, full
    from pandas import DataFrame
#    import csvreader
    
#    print('Creating network')    
    
    [A, Z] = Create_network(DF_Network, Network_headers, Wires_headers, wire_file)

    
    
    labels = [str(i+1) for i in range(len(S_n[0]))]


    df_P_PV_av = Create_PV_DF(S_n, t_0, t_simulation, P_PV_ref = P_PV_peak)
    df_P_PV = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)    

    df_P_Grid = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_P_Grid.loc[0] = df_P_L.loc[0]

    # Only for nodes with load
    df_P_BESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_P_ref = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)

    df_SOC_BESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_C_BESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)

    df_SOC_BESS = assign_DF_values(df_SOC_BESS, 0, selection, 0.5)
    df_C_BESS = assign_DF_values(df_C_BESS, 0, selection, Capacity_BESS_BoL)
    
    # For every node
    df_P_n = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_V_n = DataFrame(zeros([int(t_simulation), len(A)]))
    df_I_n = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)

    
    
    V_n_registry = empty([len(A), int(t_simulation)])
    I_registry = empty([len(A), int(t_simulation)]) 

   
    if not V_0:
        V_0 = zeros(len(A),dtype=complex) + V_feeder
        
    if not I_0:
        I_0 = zeros(len(A),dtype=complex)
#        I_0 = Estimate_Currents(A, [i[0] for i in S_n], I_0, V_0)
        I_0 = Estimate_Currents(A, [1000*S for S in df_P_Grid.loc[0].tolist()], I_0, V_0)        
            
    [V_n, I] = Network_state(A, Z, df_P_Grid, 0, I_0, V_0) # + t_0
        
    V_n_registry[:,0] = V_n
    I_registry[:,0] = I    
    
#    P_ref_min = -1.5
#    P_ref_max = 1.5
        
    print('Start Simulation')   
      

    for timestep in range(1,int(t_simulation)):
#        print('---', timestep, '---') # , 4*df_Prices['Current'][timestep]
        if (timestep)%(24/dt) == 0:
            print('Day', t_0/(24/dt) + int(timestep/(24/dt)))
        start_ts = time.time()
        

        [df_P_BESS, df_SOC_BESS, df_C_BESS, df_P_ref, df_P_Grid] = Network_Control_Centralized(A, Z, I_0, V_0, selection, remaining, df_P_L, df_P_PV_av, df_P_PV, df_P_BESS, df_SOC_BESS, df_C_BESS, df_P_ref, timestep, t_0, Capacity_BESS_BoL = 1000, V_ref_min = V_ref_min, V_ref_max = V_ref_max, P_BESS_max = P_BESS_max)
        
        
        [V_n, I] = Network_state(A, Z, df_P_Grid, timestep, I_0, V_0)
        
        V_n_registry[:,timestep] = V_n
        I_registry[:,timestep] = I
        
        
        end_ts = time.time()
        t_registry.append(end_ts-start_ts)

    df_V_n = DataFrame(V_n_registry).T
    df_V_n.columns = labels
    df_I_n = DataFrame(I_registry).T
    df_I_n.columns = labels

    return [V_n_registry, I_registry, [df_P_L, df_P_BESS, df_P_Grid, df_P_ref, df_SOC_BESS, df_C_BESS, df_P_PV_av, df_P_PV, df_P_n, df_V_n, df_I_n]]

###############################################################################

def Network_Control_Centralized(A, Z, I_0, V_0, selection, remaining, df_P_L, df_P_PV_av, df_P_PV, df_P_BESS, df_SOC_BESS, df_C_BESS, df_P_ref, timestep, t_0 = 0, Capacity_BESS_BoL = 1000, V_ref_min = 0.95*400, V_ref_max = 1.05*400, P_BESS_max = 1000):
    from pandas import DataFrame
    
    [V_n, I] = Network_state(A, Z, df_P_L, timestep, I_0, V_0)

    if (min(V_n) > V_ref_min) and (max(V_n) > V_ref_min):   # No action required
        for node in selection:        
        
            # Charge BESS
            [df_P_BESS.iloc[timestep,node], df_SOC_BESS.iloc[timestep,node], df_C_BESS.iloc[timestep,node], df_P_ref.iloc[timestep,node], df_P_PV.iloc[timestep,node]] = Central_BESS_Control(df_SOC_BESS.iloc[timestep-1,node], df_C_BESS.iloc[timestep-1,node], df_P_L.loc[timestep].sum(), df_P_PV_av.iloc[timestep,node], operation = "Charge", P_BESS_max = P_BESS_max)
        
    else:
        for node in selection:                
            #Discharge BESS
            [df_P_BESS.iloc[timestep,node], df_SOC_BESS.iloc[timestep,node], df_C_BESS.iloc[timestep,node], df_P_ref.iloc[timestep,node], df_P_PV.iloc[timestep,node]] = Central_BESS_Control(df_SOC_BESS.iloc[timestep-1,node], df_C_BESS.iloc[timestep-1,node], df_P_L.loc[timestep].sum(), df_P_PV_av.iloc[timestep,node], operation = "Support", P_BESS_max = P_BESS_max)

    
    return [df_P_BESS, df_SOC_BESS, df_C_BESS, df_P_ref, DataFrame(df_P_L.values + df_P_BESS.values, columns = df_P_L.columns)]


###############################################################################
    
def Central_BESS_Control(SoC_0, Capacity_BESS, P_Grid_sum, P_PV_av, operation = "Charge", charge_time = 15, P_BESS_max = 1000, Capacity_BESS_BoL = 1000):
    from MCES_library import BESS_perm_max, BESS_perm_min
    from numpy import clip
    
    
    if operation == "Charge":    
#        for node in selection:
        P_BESS = clip(min([-SoC_0*Capacity_BESS/charge_time, -P_PV_av]), BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max))
        [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
    
        return [-P_BESS, SoC_BESS, Capacity_BESS, 0, P_PV_av] #  clip(P_PV_av, 0, -P_BESS)
            
    elif operation == "Support":
#        for node in selection:            
#        df_S_n.loc[timestep].mean()
        P_BESS = clip(P_Grid_sum - P_PV_av, BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max))
        [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
        
        return [-P_BESS, SoC_BESS, Capacity_BESS, P_Grid_sum, P_PV_av]

###############################################################################

def Voltage_Power_correlation(all_nodes, df_P_Grid, df_V_n, df_Prices, V_feeder = 400, t0 = 0, tf = 1*24*4):
    
    from numpy import linspace, poly1d, polyfit
    import matplotlib.pyplot as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })      

    fig, axs = plt.subplots(3, 1, sharex=True) 
    
    # Power
    axs[0].grid()
    axs[0].plot(df_P_Grid[[str(i + 1) for i in all_nodes]].min(axis=1).tolist(), color='r', label='Min')
    axs[0].plot(df_P_Grid[[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), color='y', label='Mean')
    axs[0].plot(df_P_Grid[[str(i + 1) for i in all_nodes]].max(axis=1).tolist(), color='b', label='Max')
    axs[0].set_ylabel('$P_{Grid}$, [kW]')
    axs[0].set_xlim([t0,tf]) 
    axs[0].legend(loc='upper left', ncol = 1, fontsize=10)

    # Voltage
    axs[1].grid()
    axs[1].plot([j/V_feeder for j in df_V_n[[str(i + 1) for i in all_nodes]].min(axis=1).tolist()], color='r', label='Min')
    axs[1].plot([j/V_feeder for j in df_V_n[[str(i + 1) for i in all_nodes]].mean(axis=1).tolist()], color='y', label='Mean')
    axs[1].plot([j/V_feeder for j in df_V_n[[str(i + 1) for i in all_nodes]].max(axis=1).tolist()], color='b', label='Max')
    axs[1].set_ylabel('$V_{n}$, [p.u.]')
    axs[1].set_xlim([t0,tf]) 
    axs[1].legend(loc='lower left', ncol = 1, fontsize=10)


    # Price
    axs[2].grid()
    axs[2].plot([4*i for i in df_Prices['Q25'].tolist()], color='k', label='Q25')
    axs[2].plot([4*i for i in df_Prices['Q75'].tolist()], color='k', label='Q75')
    axs[2].plot([4*i for i in df_Prices['Current'].tolist()], color='b', label='Current')
    axs[2].set_xlim([t0,tf]) 
    axs[2].set_ylabel('$c_{Grid}$, [€/kWh]')
    axs[2].legend(loc='upper left', ncol = 1, fontsize=10)

    l_P = []
    l_V = []
    
    for P,V in zip(df_P_Grid[[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), [j/V_feeder for j in df_V_n[[str(i + 1) for i in all_nodes]].min(axis=1).tolist()]):
        if P > 0:        
            l_P.append(P)
            l_V.append(V)
    for P,V in zip(df_P_Grid[[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), [j/V_feeder for j in df_V_n[[str(i + 1) for i in all_nodes]].max(axis=1).tolist()]):
        if P < 0:        
            l_P.append(P)
            l_V.append(V)    
        
    model = poly1d(polyfit(l_P, l_V, 2))
    
    print('Voltage-Power correlation equation: ', model)
    
    
        
    plt.figure(constrained_layout=True)
    plt.grid()
    plt.scatter(df_P_Grid[[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), [j/V_feeder for j in df_V_n[[str(i + 1) for i in all_nodes]].min(axis=1).tolist()], color = 'r', label = 'V$_{min}$', marker = '.', alpha=0.1)
    plt.scatter(df_P_Grid[[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), [j/V_feeder for j in df_V_n[[str(i + 1) for i in all_nodes]].max(axis=1).tolist()], color = 'b', label = 'V$_{max}$', marker = '.', alpha=0.1)
    plt.plot(linspace(min(l_P), max(l_P), 130), model(linspace(min(l_P), max(l_P), 130)), color ='k', label='$\hat{V}$') # , label=eq_latex
    plt.xlabel('Average grid power, $\overline{P}_{Grid}$, [kW]')
    plt.ylabel('Voltage, V, [p.u.]')
    plt.legend(loc='lower left')
#    plt.title('Model: '+"${}$".format(eq_latex))
    plt.show()
    
    

###############################################################################

def Unify_array(arr):
    l = arr.tolist()
    l2 = []
    for i in l:
        l2 += i
        
    return l2

###############################################################################

def Make_voltage_boxplots(data, headers, V_feeder = 400): #  V_feeder = 230
    from numpy import arange, array
    import matplotlib.pylab as plt

    boxplot_ticks = [''] + headers + ['']

    plot_lists = [Unify_array(i/V_feeder) for i in data]

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })  
    
    plt.figure(constrained_layout=True)
    plt.boxplot(plot_lists)
    plt.title('Voltage distribution per method')
#    plt.xlabel('Node')
    plt.xticks(arange(0, len(data)+2, step=1), array(boxplot_ticks))      # arange(1, len(V_n)+1, step=1)
    plt.ylabel('Voltage [pu]')
    plt.grid()
    plt.show()   

###############################################################################

def count_V_compliance(Case_n, lim = 0.05, percentage = True, index = 4, V_feeder = 400):
    v_low = []
    v_up = []
    for case in Case_n:
        i_low = []
        i_up = []
        for node in case[index]:
            for ts in node:
                
                if ts<(1 - lim)*V_feeder:
                    i_low.append(ts)
                elif ts>(1 + lim)*V_feeder:
                    i_up.append(ts)
                    
        v_low.append(len(i_low))  
        v_up.append(len(i_up)) 
    
    if percentage:
        v_low = [100*i/(len(Case_n[0][index])*len(Case_n[0][index][0])) for i in v_low]
        v_up = [100*i/(len(Case_n[0][index])*len(Case_n[0][index][0])) for i in v_up]

    return [v_up, v_low]

###############################################################################

def Plot_grid_behaviour(case_n, case_centralized, node = 17, i = -1, start = 0, end = 1, penetrations = [10*i for i in arange(1,10,1)], V_feeder = 400):
    import matplotlib.pylab as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 13
    })       
        
    fig, axs = plt.subplots(3, sharex=True) # , constrained_layout=True

    ndays = (end - start)*4*24

    axs[0].grid()
    axs[0].set_ylabel('$V_{n}$, [p.u.]')
    axs[0].plot([min(case_n[i][27].loc[ts])/V_feeder for ts in range(ndays)], color='b', label='Case 2')
    axs[0].plot([max(case_n[i][27].loc[ts])/V_feeder for ts in range(ndays)], color='b')
    axs[0].plot([min(case_centralized[i][15].loc[ts])/V_feeder for ts in range(ndays)], color='r', label='Case 8')
    axs[0].plot([max(case_centralized[i][15].loc[ts])/V_feeder for ts in range(ndays)], color='r')
    axs[0].legend(loc='lower right', ncol = 1, fontsize=10)
    
    axs[1].grid()
    axs[1].set_ylabel('$P$, [kW]')
    axs[1].plot(case_centralized[i][7][str(node+1)].tolist(), color='r', label='BESS')
    axs[1].plot([-P_DSO for P_DSO in case_centralized[i][9][str(node+1)].tolist()], color='b', label='Setpoint')
    axs[1].plot(case_centralized[i][13][str(node+1)].tolist(), color='g', label='PV') 
    axs[1].legend(loc='lower right', ncol = 3, fontsize=10)

    
    
    axs[2].grid()
    axs[2].set_ylabel('$SoC$, [%]')
    axs[2].plot([100*SoC for SoC in case_centralized[i][10][str(node+1)].tolist()], color='r')
    axs[2].set_xlim([0,ndays])

###############################################################################

def case_parameters(case):
    if case == 1 or case == 0:
        enable_HP = False
        BESS_power = 0      
        BESS_energy = 10
        enable_TESS = False
        TESS_capacity = 4
        follow_control = False
    elif case == 2:
        enable_HP = True
        BESS_power = 0      
        BESS_energy = 10
        enable_TESS = False
        TESS_capacity = 4     
        follow_control = False
    elif case == 3:
        enable_HP = True
        BESS_power = 10     
        BESS_energy = 10
        enable_TESS = False
        TESS_capacity = 4 
        follow_control = False
    elif case == 4:
        enable_HP = True
        BESS_power = 10     
        BESS_energy = 10
        enable_TESS = True
        TESS_capacity = 4
        follow_control = False
    elif case == 5:
        enable_HP = True
        BESS_power = 10     
        BESS_energy = 10
        enable_TESS = True
        TESS_capacity = 4
        follow_control = True   
    elif case == 6:
        enable_HP = True
        BESS_power = 10      
        BESS_energy = 10
        enable_TESS = True
        TESS_capacity = 4
        follow_control = False        
    elif case == 7:
        enable_HP = True
        BESS_power = 10      
        BESS_energy = 10
        enable_TESS = False
        TESS_capacity = 4 
        follow_control = True    
        

    return [enable_HP, BESS_power, BESS_energy, enable_TESS, TESS_capacity, follow_control]  
     
###############################################################################

def TESS_SoC(T_TESS, T_min = 50+273, T_max = 90+273):
    
    return (100/(T_max - T_min))*T_TESS + (-(100/(T_max - T_min))*T_min)

###############################################################################
            
def Plot_house_behaviour(case_n, node, penetration = 20, case = 0, i = -1, start = 0, end = 24*4*7, penetrations = [10*i for i in arange(1,10,1)]):
    from numpy import array#linspace, sin, pi
    import matplotlib.pylab as plt

    plt.rcParams.update({
        "font.family": "Times New Roman",
        'font.size': 13
    })      

    if i != -1:
        if penetration != 100:
            i = 5*(penetrations.index(penetration)) + case
        else:
            i = -1
    
    P_L = case_n[i][8][str(node)].tolist()
    P_PV_av = case_n[i][6][str(node)].tolist()
    P_PV = case_n[i][7][str(node)].tolist()
    P_HP = case_n[i][9][str(node)].tolist()
    P_HP_TESS = case_n[i][10][str(node)].tolist()
    P_BESS = case_n[i][11][str(node)].tolist()
    SoC_BESS = [100*i for i in case_n[i][15][str(node)].tolist()]
    P_Grid = case_n[i][12][str(node)].tolist()
    P_Grid_DSO = case_n[i][13][str(node)].tolist()
    Prices = case_n[i][29]['Current'].tolist()
    
    Qdot_L = [i/1000 for i in case_n[i][17][str(node)].tolist()]
    Qdot_Boiler = [i/1000 for i in case_n[i][23][str(node)].tolist()]
    Qdot_HP = [i/1000 for i in case_n[i][18][str(node)].tolist()]
    Qdot_HP_TESS = [i/1000 for i in case_n[i][19][str(node)].tolist()]
    Qdot_TESS = [i/1000 for i in case_n[i][20][str(node)].tolist()]
    
    T_in = case_n[i][25][str(node)].tolist()
    T_set_day = [273 + 20 - 3]*int((6-0)*4) + [273 + 20]*int((22-6)*4) + [273 + 20 - 3]*int((24-22)*4)
    T_set = array(T_set_day*int(len(T_in)/24/4)) 
    SoC_TESS = [TESS_SoC(i) for i in case_n[i][21][str(node)].tolist()]
    
    
    fig, axs = plt.subplots(7, 2, sharex=True) # , constrained_layout=True

    # Electric
    axs[0, 0].grid()
    axs[0, 0].plot(P_L) # Load
    axs[0, 0].set_ylabel('$P_{L}$, [kW]')
    axs[1, 0].grid()
    axs[1, 0].plot(P_PV_av, color='r', label='$P_{av}$') # PV av
    axs[1, 0].plot(P_PV, color='b', label='$P_{PV}$') # PV av, PV
    axs[1, 0].set_ylabel('$P_{PV}$, [kW]')
    axs[1, 0].legend(loc='upper left', ncol = 2, fontsize=10)
    axs[2, 0].grid()
    axs[2, 0].plot(P_HP_TESS, color='b', label='TESS') # TESS    
    axs[2, 0].plot(P_HP, color='r', label='Load') # HP    
    axs[2, 0].set_ylabel('$P_{HP}$, [kW]')
    axs[2, 0].legend(loc='upper right', ncol = 1, fontsize=10)
    axs[3, 0].grid()
    axs[3, 0].plot(P_BESS) # BESS
    axs[3, 0].set_ylabel('$P_{BESS}$, [kW]')
    axs[4, 0].grid()
    axs[4, 0].plot(SoC_BESS) # BESS
    axs[4, 0].set_ylabel('$SoC_{BESS}$, [%]')
    axs[4, 0].set_ylim([15,95])    
    axs[5, 0].grid()
    axs[5, 0].plot(P_Grid, label = 'Grid')
    axs[5, 0].plot(P_Grid_DSO, label = 'DSO')
    axs[5, 0].set_ylabel('$P_{Grid}$, [kW]')
    axs[5, 0].legend(loc='upper left', ncol = 2, fontsize=10)
    axs[6, 0].grid()
    axs[6, 0].plot(Prices)
    axs[6, 0].set_ylabel('$\lambda$, [€/kWh]')
    axs[6, 0].set_xlim([start,end])
    
    # Thermal
    axs[0, 1].grid()
    axs[0, 1].plot(Qdot_L)
    axs[0, 1].set_ylabel('$\dot{Q}_{L}$, [kW]')    
    axs[1, 1].grid()
    axs[1, 1].plot(Qdot_Boiler)
    axs[1, 1].set_ylabel('$\dot{Q}_{B}$, [kW]')    
    axs[2, 1].grid()    
    axs[2, 1].plot(Qdot_HP_TESS, color='b', label='TESS') # TESS    
    axs[2, 1].plot(Qdot_HP, color='r', label='Load') # HP
    axs[2, 1].set_ylabel('$\dot{Q}_{HP}$, [kW]')  
    axs[2, 1].legend(loc='upper right', ncol = 1, fontsize=10)    
    axs[3, 1].grid()
    axs[3, 1].plot(Qdot_TESS)
    axs[3, 1].set_ylabel('$\dot{Q}_{TESS}$, [kW]') 
    axs[4, 1].grid()
    axs[4, 1].plot(SoC_TESS) # BESS
    axs[4, 1].set_ylabel('$SoC_{TESS}$, [%]')
    axs[4, 1].set_ylim([0,100])    
    axs[5, 1].grid()
    axs[5, 1].plot([i-273 for i in T_in], color='r', label='in')
    axs[5, 1].plot([i-273 for i in T_set], color='b', label='set')
    axs[5, 1].set_ylabel('T, [°C]')
    axs[5, 1].legend(loc='lower right', ncol = 1, fontsize=10)      
    axs[6, 1].grid()
    axs[6, 1].plot(Prices)    
    axs[6, 1].set_ylabel('$\lambda$, [€/kWh]')
    axs[6, 1].set_xlim([start,end])


###############################################################################
###########################   SciPy Optimization   ############################


# Define the objective function to minimize (profit function)
def Cost_function_1(x, *args): #, *args
    # from numpy import divide, matmul, transpose
    Sn0 = args[5]    
    return sum([(P - Sn)**2 for P, Sn in zip(x, Sn0)])  

###############################################################################
    
def Cost_function_2(x, *args): #, *args
    from numpy import divide, matmul, transpose
    A = args[0]
    B = args[1]
    Z = args[2]
    V0 = args[3]
    I0 = args[4]  
    return sum([abs(abs(matmul(A[i,:],V0) - matmul(Z[i,:],(matmul(transpose(A),transpose(I0)) + divide(x, [i/1000 for i in V0]))) + B[i]*400) - 400)**2 for i in range(301)])

###############################################################################

def Cost_function(x, *args):
    
    weights = args[6]
    
    return weights[0]*Cost_function_1(x, *args) + weights[1]*Cost_function_2(x, *args)

###############################################################################

def Find_power_setpoints(A, Z, B, V0, I0, Sn0, all_nodes, selected_nodes, pu_lim = 0.05, v_feeder = 400, weights = [1162690, 13], plotting = False):
    from scipy.optimize import minimize
    # https://www.youtube.com/watch?v=X0LvnxSqfNk&ab_channel=KodyPowell
    
    # Define the bounds for each decision variable
    
    # Initial guess for the decision variables
    initial_guess = [i for i in Sn0]
    
    # Define the inequality constraint functions    

    cons = []
    c1_low = {'type': 'ineq', 'fun': lambda x : -sum(x) + 161}
    cons.append(c1_low)
    c1_high = {'type': 'ineq', 'fun': lambda x : sum(x) + 175} # 1.05 -> 175m 1.03 -> 105
    cons.append(c1_high)    
    for factor in range(301):   # Check for the nodes without load
        # Define the equality constraint functions
        # Distribution nodes
        if factor not in all_nodes:
            c2 = {'type': 'eq', 'fun': lambda x, i=factor : x[i] - 0}
            cons.append(c2)
        # Consumer nodes
        elif factor not in selected_nodes:
            c3 = {'type': 'eq', 'fun': lambda x, i=factor : x[i] - Sn0[i]}
            cons.append(c3)            
    
    opt = minimize(Cost_function_1, initial_guess, method = 'SLSQP', args = (A, B, Z, V0, I0, Sn0, weights), constraints = cons, options={'disp': True}) #, method = 'SLSQP', bounds = bounds, constraints = cons, args = (A, Z, df_S_n_0, ts, I0, V0), options={'maxiter': 1000}

    if plotting:
        import matplotlib.pyplot as plt 
        [V_n, I_0] = Estimate_Node_voltage(A, Z, [i*1000 for i in opt.x], I0, V0, show_plot = True)       
                
        fig, axs = plt.subplots(2, sharex=True)
    # Power    
        axs[0].plot(Sn0, color = 'b', label = 'Sn0')
        axs[0].plot([i for i in opt.x], color = 'r', label = 'Opt.x')
        axs[0].set_ylabel('Node power, $S_{n}$, [kW]')
        axs[0].legend(loc='upper right')
        axs[0].grid()    

    # Voltage
        axs[1].plot([i/400 for i in V0], color = 'b', label = 'V0')
        axs[1].plot([i/400 for i in V_n], color = 'r', label = 'Opt.x')
        axs[1].set_ylabel('Node voltage, $V_{n}$, [p.u.]')
        axs[1].legend(loc='upper right', ncol = 2)
        axs[1].set_xlim([0,301])
        axs[1].grid() 
        axs[1].set_xlabel('Node')     


        
    return opt.x