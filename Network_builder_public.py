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

###############################################################################


def Clear_zeros_row_col(arr):
    from numpy import argwhere, all, delete
    
    arr = arr[~all(arr == 0, axis=1)]    
    idx = argwhere(all(arr[..., :] == 0, axis=0))    
    return delete(arr, idx, axis=1)
#    return arr

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
    
#        S_n*=1000
        if show_plots:
            import matplotlib.pylab as plt
            plt.rcParams.update({
            #    "text.usetex": True,
                "font.family": "Times New Roman",
                'font.size': 16
            })
    
    #        labels = [str(i[1]) for i in Load_nodes]
            
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
    
#    Load_nodes = []
#    
#    for i, node, profile in zip(DF_Network[Network_headers[4]].isnull(), DF_Network[Network_headers[1]], DF_Network[Network_headers[4]]):
#        if not i:
#            Load_nodes.append([node, int(findall(r"\d+", profile)[0])])
#    
#    
#    for i in range(len(Load_nodes)):
#        Load_nodes[i].append(Create_profile(True, Load_nodes[i][1]))
#    
#    if show_plots:
#        import matplotlib.pylab as plt
#        plt.rcParams.update({
#        #    "text.usetex": True,
#            "font.family": "Times New Roman",
#            'font.size': 16
#        })
#
##        labels = [str(i[1]) for i in Load_nodes]
#        
#        plt.figure(constrained_layout=True)  
#        plt.grid()
#        for i in Load_nodes:
#            plt.plot(i[2], label = 'Node '+str(i[0])+', '+str(i[1])+' kWh')
#        plt.ylabel('Demanded power, $P_{Grid}$, [kW]')
#        plt.xlabel('Time, t, [s]') 
#        plt.legend()
#        plt.show()
#  
#
#    return Load_nodes

###############################################################################   

def Create_admittance_matrix(DF_Network, Network_headers, show_heatmap = False):
    from numpy import zeros
    
    A = zeros((max(DF_Network[Network_headers[1]]), max(DF_Network[Network_headers[1]]))) 
    
    for i,j in zip(DF_Network[Network_headers[0]],DF_Network[Network_headers[1]]):
#        A[i-1,j-1] = 1
        A[j-1,i-1] = 1
        
#    A = Clear_zeros_row_col(A)
    
    if show_heatmap: 
#        from numpy import arange, all
        import matplotlib.pylab as plt
        plt.rcParams.update({
        #    "text.usetex": True,
            "font.family": "Times New Roman",
            'font.size': 16
        })

        
        plt.figure(constrained_layout=True)
#        fig, (ax) = plt.subplots(1, 1, constrained_layout=True, sharey=True)
#        a = plt.imshow(A, cmap = 'binary')
#        cbar = fig.colorbar(a, ax=ax)
#        cbar.set_label('Admittance, [-]')
        plt.imshow(A, cmap = 'binary')
#        plt.xticks(arange(0, len(A), step=1), arange(1, len(A)+1, step=1))
#        plt.yticks(arange(0, len(A), step=1), arange(1, len(A)+1, step=1))
        plt.title('Admittance Matrix')
        plt.xlabel('Node')
        plt.ylabel('Node')
        plt.show()        
    
    return A#[~all(A == 0, axis=1)]

###############################################################################

def Create_adjacency_matrix(selected_nodes):
    from numpy import zeros

    A = zeros((len(selected_nodes), len(selected_nodes)))   
    
    for i in range(len(selected_nodes)-1):
        A[i+1,i] = 1
        A[i,i+1] = 1
    
    return A

###############################################################################

def Create_degree_matrix(A):
    from numpy import zeros
    
    D = zeros((len(A), len(A)))    
    for d in range(len(A)):
        D[d,d] = sum(A[d,:]) + sum(A[:,d])
        
    return D

###############################################################################

def Create_Laplacian_matrix(A, D):
    return D - A

###############################################################################
    
def Calculate_concensus(A, D, x, step_size = 0.25, max_steps = 100, plotting = False):
    from numpy import array, zeros
    
    x=array(x)
    
    x_new = zeros(len(x))
    x_registry = zeros((len(x), max_steps))
    
    L = Create_Laplacian_matrix(A, D)
    
    for step in range(max_steps):
        x_registry[:,step] = x
#        for i in range(len(x)):
#            x_new[i] = x[i] + step_size*sum([A[j,i]*(x[j] - x[i]) for j in range(len(x))])
#        
#        x = [i for i in x_new]
        x_new = x - step_size*L.dot(x)
        x = x_new
    
    
    
    if plotting:
        import matplotlib.pylab as plt
        plt.rcParams.update({
        #    "text.usetex": True,
            "font.family": "Times New Roman",
            'font.size': 16
        })
        plt.figure(constrained_layout=True)
        for x in x_registry:
            plt.plot(x)
        
    return x_registry
    
###############################################################################
    
def Concensus_Example():
    from numpy import ones

    A = ones((5,5))
#    D = zeros((5,5))
    
    for i in range(len(A)):
        A[i,i] = 0
    
    x_0 = [4,3,1,3,1]
    
    return Calculate_concensus(A, 0, x_0, plotting = True) 
    

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
#        Z[i-1,j-1] = Calculate_impedance(DF_Wires, k, l, Wires_headers, km)
        Z[j-1,i-1] = Calculate_impedance(DF_Wires, k, l, Wires_headers, km)

#    Z = Clear_zeros_row_col(Z)
    if show_heatmap:
#        import seaborn as sns
#        import matplotlib.pylab as plt
#        
#        
#        plt.figure(constrained_layout=True)
#        sns.heatmap(abs(Z))
#        plt.title('Impedance Matrix')
#        plt.xlabel('Node')
#        plt.ylabel('Node')
#        plt.show()       
#        from numpy import arange        
        import matplotlib.pylab as plt
        plt.rcParams.update({
        #    "text.usetex": True,
            "font.family": "Times New Roman",
            'font.size': 16
        }) 

        
#        plt.figure(constrained_layout=True)
#        plt.imshow(abs(Z), cmap = 'binary')
        
        fig, (ax) = plt.subplots(1, 1, sharey=True) # , constrained_layout=True
        a = plt.imshow(abs(Z), cmap = 'binary')
        cbar = fig.colorbar(a, ax=ax)
        cbar.set_label('Impedance, [$\Omega$]')        
        
#        plt.xticks(arange(0, len(Z), step=1), arange(1, len(Z)+1, step=1))
#        plt.yticks(arange(0, len(Z), step=1), arange(1, len(Z)+1, step=1))
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
    #    "text.usetex": True,
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
    from numpy import array, divide, matmul, transpose, where, isnan
    
#    from pandas import isna

#    print('----------------------------------')
#    print(type(A))
#    print(type(S_n))    
    
    if [] != where(isnan(divide(S_n, V_0))):
#    if [] != where(isna(divide(S_n, V_0))):        
        I_n = divide(S_n, V_0)
    else:
        print('zeros')
        I_n = array([S/V if V!=0 else 0 for S,V in zip(S_n, V_0)])

#    I_n = divide(S_n, V_0)

    I = matmul(transpose(A),transpose(I_0)) + I_n
#    I = [a+b for a,b in zip(matmul(transpose(A),transpose(I_0)), I_n)]
#    I = I.tolist()[1:]
    
#    print('----------------------------------')
#    print(I)
#    print(type(I_0))
#    print(I_0)
    
    i = 0
    while max(abs(I - I_0)) > dI:
        
        I_0 = I
        I = matmul(transpose(A),I_0) + I_n
        
#        i+=1
        if i>100:
            print('Error in current estimation')
            break

    if show_plot:
        plot_current(I)
        
    return I

###############################################################################

def plot_voltage(V_n, V_feeder = 400): # V_feeder = 230
    from numpy import arange
    import matplotlib.pylab as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })  
    
    plt.figure(constrained_layout=True)
    plt.plot([abs(i)/V_feeder for i in V_n])
    plt.title('Voltage vector')
    plt.xlabel('Node')
    plt.xlim([1, len(V_n)-1])
#    plt.xticks(arange(0, len(V_n), step=1), arange(1, len(V_n)+1, step=1))      
    plt.ylabel('Voltage [pu]')
    plt.grid()
    plt.show()   

###############################################################################
    
def plot_power(S_n):
    from numpy import arange
    import matplotlib.pylab as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })  
    
    plt.figure(constrained_layout=True)
#    plt.plot([abs(i) for i in S_n])
    plt.plot([i/1000 for i in S_n])    
    plt.title('Power vector')
    plt.xlabel('Node')
    plt.xlim([1, len(S_n)-1])
    plt.xticks(arange(0, len(S_n), step=1), arange(1, len(S_n)+1, step=1))      
    plt.ylabel('Power [kVA]')
    plt.grid()
    plt.show()       

###############################################################################

def Estimate_Node_voltage(A, Z, S_n, I_0, V_0, V_feeder = 400, dV = 0.001, dI = 0.01, show_plot = False): # V_feeder = 400
    from numpy import matmul, zeros
    
    I_0 = Estimate_Currents(A, S_n, I_0, V_0, dI)
    
    B = zeros(len(I_0))
    B[0] = 1

    V_n = matmul(A,V_0) - matmul(Z,I_0) + B*V_feeder
    
#    V_n = zeros(len(I_0)) + V_feeder
    
#    print(I_0)
#    print(V_n)
    
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
     
#    print(I_0)
        
    return [V_n, I_0]

###############################################################################

def Create_network(DF_Network, Network_headers, Wires_headers, wire_file = 'Gaia_cables.xlsx'):
    A = Create_admittance_matrix(DF_Network, Network_headers, False)
    Z = Create_impedance_matrix(DF_Network, Network_headers, Wires_headers, False, wire_file)

    return [A, Z]    
    
###############################################################################

def Network_state(A, Z, df_S_n, timestep, I_0, V_0):
    
#    S_n_timestep = S_n[timestep,:]
    S_n_timestep = df_S_n.loc[timestep].tolist()
    return Estimate_Node_voltage(A, Z, [1000*S for S in S_n_timestep], I_0, V_0)
    

###############################################################################
    
def Create_PV_DF(S_n = False, t_0 = 0, t_simulation = 24*4, H = 0, P_PV_ref = False, module_power_ref = 0.315, return_modules = False, dt = 0.25):
    import csvreader
    from pandas import DataFrame
    
    CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
    CSVDataPV.data2array()
    labels = [str(i+1) for i in range(len(S_n[0]))]
    

    if P_PV_ref != False:#P_PV_ref.any():
        P_PV = DataFrame([[ref*i[0]/module_power_ref/1000 for i in CSVDataPV.ar[t_0:t_0+t_simulation+H]] for ref in P_PV_ref])    

    elif P_PV_ref == False:#S_n.any():
        E_PV_ref_module = sum([i[0] for i in CSVDataPV.ar])*dt/1000
#        print([sum(S_n[:,node]) for node in range(len(S_n[0]))])
        
        n_modules = [abs(sum(S_n[:,node])*dt//E_PV_ref_module) for node in range(len(S_n[0]))]
        
#        print(n_modules)
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
    
#global operation_registry 
#operation_registry = [0]

def Electric_Power_Balance(P_Load, P_HP, P_BESS, P_BESS_max, P_BESS_min, P_PV, P_grid_DSO, top = True, exact_match = True):
#    from numpy import clip
    
#    print(P_grid_DSO)
    
    if P_grid_DSO:
        if exact_match:
            if P_HP + P_Load - P_PV - P_BESS_max >= P_grid_DSO:         # The BESS and the PV cannot compensate (increase injection)
#                print('Exact match 1')
                return [P_PV, P_BESS_max, P_HP + P_Load - P_PV - P_BESS_max]

            elif (P_HP + P_Load - P_PV - P_BESS_max <= P_grid_DSO) and (P_HP + P_Load - P_PV >= P_grid_DSO):       # Curtailing the PV is not enough to compensate, but the BESS can compensate (increase injection)
#                print('Exact match 2')
                return [P_PV, P_HP + P_Load - P_PV - P_grid_DSO, P_grid_DSO]

            elif (P_HP + P_Load - P_PV <= P_grid_DSO) and (P_HP + P_Load >= P_grid_DSO):                    # Curtailing the PV is enough to compensate
#                print('Exact match 3')
                return [P_HP + P_Load - P_grid_DSO, 0, P_grid_DSO]

            elif P_HP + P_Load - P_BESS_min >= P_grid_DSO:              # The BESS has enough space to consume power to compensate (increase consumption)
#                print('Exact match 4')
                return [0, P_HP + P_Load - P_grid_DSO, P_grid_DSO]

            elif P_HP + P_Load - P_BESS_min <= P_grid_DSO:               # The BESS cannot consume enough power to compensate (increase consumption)          
#                print('Exact match 5')
                return [0, P_BESS_min, P_HP + P_Load - P_BESS_min]        
            
           
        else:
            if top:     # Peak-Shaving         
                if P_HP + P_Load - P_PV <= P_grid_DSO:      # No action required
    #                operation_registry.append(1)
                    return [P_PV, 0,  P_Load + P_HP - 0 - P_PV] 
                
                elif P_HP + P_Load - P_PV - P_BESS_max <= P_grid_DSO:   # BESS can compensate it
    #                operation_registry.append(2)
                    return [P_PV, P_HP + P_Load - P_PV - P_grid_DSO, P_grid_DSO]
                
                else:    # BESS cannot compensate completely
    #                operation_registry.append(3)
                    return [P_PV, P_BESS_max, P_Load + P_HP - P_BESS_max - P_PV]
        
            else:       # Power curtailment     
                if P_HP + P_Load - P_PV >= P_grid_DSO:     # No action required
    #                operation_registry.append(4)
                    return [P_PV, 0, P_Load + P_HP - 0 - P_PV]
                
                elif P_HP + P_Load >= P_grid_DSO:   # PV curtailment is enough
    #                operation_registry.append(5)
                    return [P_HP + P_Load - P_grid_DSO, 0, P_grid_DSO]
                
                elif P_HP + P_Load + P_BESS_min >= P_grid_DSO:    # BESS can compensate, curtailment is not enough
    #                operation_registry.append(6)
                    return [0, P_grid_DSO - P_HP - P_Load, P_grid_DSO]            
                
                else:   # BESS cannot compensate
    #                operation_registry.append(7)
                    return [0, P_BESS_min, P_Load + P_HP - P_BESS_min - P_PV]            
                
    else:           # No request from the DSO   
#        operation_registry.append(8)
        return [P_PV, P_BESS, P_Load + P_HP - 0 - P_PV]     


###############################################################################

#global control_registry 
#control_registry = [0]

def Check_control_cases(T_in, T_set, T_amb, T_TESS_0, T_ret, T_soil, G, P_Load, P_grid_DSO, top, SoC_0, Cost_grid, Cost_grid_median, Cost_grid_Q1, Cost_grid_Q3, Capacity_BESS, P_PV_av, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP = True, enable_TESS = True, enable_HP_2_TESS = True, external_control = False, SoC_BESS_min = 0.9, SoC_BESS_max = 0.2, m = 4000, TESS_min_op_T = 55 + 273, TESS_max_op_T = 90 + 273, T_TESS_min = 50 + 273, T_TESS_max = 90 + 273, m_dot = 0.22, c_f = 4200, dt = 0.25, House_type = 'apartment', return_values = False):
    from tabulate import tabulate
    
    enable_HP_2_TESS = enable_TESS
    
    [P_PV, P_BESS, P_HP, T_TESS, P_Grid, SoC_BESS, Capacity_BESS, Q_D, Qdot_HP, Qdot_TESS, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, Qdot_Boiler, enable_HP_2_TESS, T_in] = heuristic_control(T_in, T_set, T_amb, T_TESS_0, T_ret, T_soil, G, P_Load, P_grid_DSO, top, SoC_0, Cost_grid, Cost_grid_median, Cost_grid_Q1, Cost_grid_Q3, Capacity_BESS, P_PV_av, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP = enable_HP, enable_TESS = enable_TESS, enable_HP_2_TESS = enable_HP_2_TESS, external_control = external_control, SoC_BESS_min = 0.9, SoC_BESS_max = 0.2, m = 4000, TESS_min_op_T = 55 + 273, TESS_max_op_T = 90 + 273, T_TESS_min = 50 + 273, T_TESS_max = 90 + 273, m_dot = 0.22, c_f = 4200, dt = 0.25, House_type = House_type)

#    data_in=[
#
#    ['T_in','= ', T_in-273, ' °C'],
#    ['T_set','= ', T_set-273, ' °C'],
#    ['T_amb','= ', T_amb-273, ' °C'],
#    ['T_TESS_0','= ', T_TESS_0-273, ' °C'],
##    ['T_ret','= ', T_ret-273, ' °C'],
##    ['T_soil','= ', T_soil[0], ' °C'],
##    ['G','= ', G, ' kW/m2'],
#    ['P_Load','= ', P_Load, ' kW'],
#    ['P_grid_DSO','= ', P_grid_DSO, ' kW'],
#    ['top','= ', top, ''],
#    ['SoC_0','= ', SoC_0, ' %'],
#    ['Cost_grid','= ', Cost_grid, ' €/kWh'],
##    ['Cost_grid_median','= ', Cost_grid_median, ' €/kWh'],
#    ['Cost_grid_Q1','= ', Cost_grid_Q1, ' €/kWh'],
#    ['Cost_grid_Q3','= ', Cost_grid_Q3, ' €/kWh'],
##    ['Capacity_BESS','= ', Capacity_BESS, ' kWh'],
#    ['P_PV_av','= ', P_PV_av, ' kW'],
#    ['enable_HP','= ', enable_HP, ''],
#    ['enable_TESS','= ', enable_TESS, ''],
#    ['enable_HP_2_TESS','= ', enable_HP_2_TESS, ''],
#    ['external_control','= ', external_control, ''],
##    ['House_type','= ', House_type, ''],
#    ]        
#    data_out=[    
#    ['P_PV','= ', P_PV, ' kW'],
#    ['P_BESS','= ', P_BESS, ' kW'],
#    ['P_HP','= ', P_HP, ' kW'],
#    ['T_TESS','= ', T_TESS-273, ' °C'],
#    ['P_Grid','= ', P_Grid, ' kW'],
#    ['SoC_BESS','= ', SoC_BESS, ' '],
##    ['Capacity_BESS','= ', Capacity_BESS, ' kWh'],
#    ['Q_D','= ', Q_D/1000, ' kW'],
#    ['Qdot_HP','= ', Qdot_HP/1000, ' kW'],
#    ['Qdot_TESS','= ', Qdot_TESS/1000, ' kW'],
#    ['T_TESS_new','= ', T_TESS_new-273, ' °C'],
#    ['Qdot_TESS','= ', Qdot_TESS/1000, ' kW'],
#    ['Qdot_TESS_SD','= ', Qdot_TESS_SD/1000, ' kW'],
#    ['Qdot_Boiler','= ', Qdot_Boiler/1000, ' kW'],
#    ['enable_HP_2_TESS','= ', enable_HP_2_TESS, ' '],
#    ['T_in','= ', T_in, ' °C']]
    
    if top:
        top = 'Peak-shaving'
    else:
        top = 'Curtailment'
    
    data = [
    ['Cost_grid','= ', Cost_grid, ' €/kWh'],            
    ['Cost_grid_Q1','= ', Cost_grid_Q1, ' €/kWh'],
    ['Cost_grid_Q3','= ', Cost_grid_Q3, ' €/kWh'],
    ['', '', '', ''], 
          
    ['External_control','= ', external_control, ''],
    ['Mode','= ', top, ''],
    ['enable_HP','= ', enable_HP, ''],
    ['enable_TESS','= ', enable_TESS, ''],
#    ['enable_HP_2_TESS','= ', enable_HP_2_TESS, ''],
    ['', '', '', ''],   
      
    
    ['T_in','= ', T_in-273, ' °C'],
    ['T_set','= ', T_set-273, ' °C'],
    ['', '', '', ''],   
    

    ['P_grid_DSO','= ', P_grid_DSO, ' kW'],
    ['P_Grid','= ', P_Grid, ' kW'],
    ['P_Load','= ', P_Load, ' kW'],
    ['P_PV_av','= ', -P_PV_av, ' kW'],
    ['P_PV','= ', -P_PV, ' kW'],
    ['P_BESS','= ', P_BESS, ' kW'],
    ['P_HP','= ', P_HP, ' kW'],
    ['', '', '', ''],   
    ['SoC_0','= ', SoC_0, ' %'],
    ['SoC_BESS','= ', SoC_BESS, ' '],
    ['', '', '', ''],



    ['T_amb','= ', T_amb-273, ' °C'],
    ['Q_D','= ', Q_D/1000, ' kW'],
    ['Qdot_HP','= ', Qdot_HP/1000, ' kW'],
    ['Qdot_TESS','= ', Qdot_TESS/1000, ' kW'],
    ['Qdot_Boiler','= ', Qdot_Boiler/1000, ' kW'],
    ['', '', '', ''],   
    ['T_TESS_0','= ', T_TESS_0-273, ' °C'],   
    ['T_TESS_new','= ', T_TESS_new-273, ' °C'],
    ['Qdot_TESS_SD','= ', Qdot_TESS_SD/1000, ' kW']
    ]    


    # Creating the table
#    table_in = tabulate(data_in, headers = ['Parameter', '', 'Value', 'Unit']) # , colalign=('left', 'center', 'center', 'center')
#    data_out = tabulate(data_in, headers = ['Parameter', '', 'Value', 'Unit']) # , colalign=('left', 'center', 'center', 'center')
    table = tabulate(data, headers = ['Parameter', '', 'Value', 'Unit']) # , colalign=('left', 'center', 'center', 'center')
    
    # Printing the table
#    print(table_in)
#    print('')
#    print(data_out)
#    print('')
    print(table)    

#    if return_values:
#        return [i[2] for i in data]
    
#heuristic_control
#Check_control_cases(T_in = 16 + 273,
#                 T_set = 17 + 273,
#                 T_amb = 10 + 273,
#                 T_TESS_0 = 60 + 273,
#                 T_ret = 38 + 273,
#                 T_soil = T_soil[0],
#                 G = 0,
#                 P_Load = 1.5,
#                 P_grid_DSO = 1,
#                 top = False,
#                 SoC_0 = 0.50,
#                 Cost_grid = 0.25,
#                 Cost_grid_median = 0.15,
#                 Cost_grid_Q1 = 0.10,
#                 Cost_grid_Q3 = 0.20,
#                 Capacity_BESS = 10,
#                 P_PV_av = 0.2,
#                 enable_HP = True,
#                 enable_TESS = True,
##                 enable_HP_2_TESS = True,
#                 external_control = True,
#                 House_type = 'apartment')


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

#    [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb, T_ret)[0]
    #P_PV = curtailment*PV_out(T_amb, G*n_PV_modules)
#    P_BESS = BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)
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
#        print('T_min')
        return [False, True, False]

###############################################################################
                        
def heuristic_control(T_in, T_set, T_amb, T_TESS_0, T_ret, T_soil, G, P_Load, P_grid_DSO, top, SoC_0, Cost_grid, Cost_grid_median, Cost_grid_Q1, Cost_grid_Q3, Capacity_BESS, P_PV_av, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP = True, enable_TESS = True, enable_HP_2_TESS = True, external_control = False, SoC_BESS_min = 0.9, SoC_BESS_max = 0.2, m = 4000, TESS_min_op_T = 55 + 273, TESS_max_op_T = 90 + 273, T_TESS_min = 50 + 273, T_TESS_max = 90 + 273, T_in_crit = 16 + 273, m_dot = 0.22, c_f = 4200, dt = 0.25, House_type = 'apartment', force_temperature = False, glazing_type = 'double', Cavity_Filling = False, Cavity_type = False, Wall_Cover = False, LRoof = 0.2, LCavity_wall = 0.05):
    from MCES_library import BESS_perm_max, BESS_perm_min, HP_Power, new_house_Temperature, House_Thermal_Losses # , update_TESS
#    from EMS_comparison import Electric_Power_Balance, update_BESS_aeging
    
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
#                [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb, T_ret)
#                #P_PV = curtailment*PV_out(T_amb, G*n_PV_modules)
#                P_BESS = BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)
#                [P_PV, P_BESS, P_Grid] = Electric_Power_Balance(P_Load, P_HP + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_BESS, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), 0, False, top)
#                [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1])                               
#                [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
               # print('No external control - Heating required - Negative prices, charge BESS')
    #            control_registry.append(1)
    
            elif Cost_grid <= Cost_grid_Q1:     # Low prices, charge BESS
                TESS_active = enable_TESS*False
                HP_active = enable_HP*True
                HP_charge_TESS = enable_TESS*enable_HP*(not HP_active)*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4, Capacity_BESS_BoL, P_BESS_max, T_TESS_max)
#                [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb, T_ret)
#                #P_PV = curtailment*PV_out(T_amb, G*n_PV_modules)
#                P_BESS = BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4
#                [P_PV, P_BESS, P_Grid] = Electric_Power_Balance(P_Load, P_HP + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_BESS, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), P_PV_av, False, top)
#                [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1])                               
#                [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
               # print('No external control - Heating required - Low prices, charge BESS')
    #            control_registry.append(2)
            
    
            elif Cost_grid <= Cost_grid_Q3:     # Average prices, BESS compensate
                TESS_active = enable_TESS*True
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*False*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_BESS_ref, P_PV, setpoint] = check_P_BESS_setpoint(SoC_0, Capacity_BESS, P_Load, HP_Power(HP_active, T_amb, T_ret)[0], HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_PV_av, P_grid_DSO, P_BESS_max = P_BESS_max)
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, P_BESS_ref, Capacity_BESS_BoL, P_BESS_max, T_TESS_max, setpoint = True)        

#                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, 0, Capacity_BESS_BoL, P_BESS_max, T_TESS_max)                

                # --- #
#                [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb, T_ret)
#                #P_PV = curtailment*PV_out(T_amb, G*n_PV_modules)
#                P_BESS = 0
#                [P_PV, P_BESS, P_Grid] = Electric_Power_Balance(P_Load, P_HP + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_BESS, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), P_PV_av, P_grid_DSO, top)
#                [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1])                               
               # [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)            
               # print('No external control - Heating required - Average prices, BESS compensate')
    #            control_registry.append(3)
                
            else:                                   # High prices, discharge BESS
                TESS_active = enable_TESS*True
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*False*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4, Capacity_BESS_BoL, P_BESS_max, T_TESS_max)                
#                [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb, T_ret)
#                #P_PV = curtailment*PV_out(T_amb, G*n_PV_modules)
#                P_BESS = BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4
#                [P_PV, P_BESS, P_Grid] = Electric_Power_Balance(P_Load, P_HP + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_BESS, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), P_PV_av, False, top)
#                [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1])                               
#                [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
               # print('No external control - Heating required - High prices, discharge BESS')
    #            control_registry.append(4)
                
            
        else:               # Heating not required
            if (Cost_grid <= 0)*(not external_control): # Negative prices, charge BESS
                TESS_active = enable_TESS*False
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*(T_TESS_0 < T_TESS_max)*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, 0, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), Capacity_BESS_BoL, P_BESS_max, T_TESS_max)                
#                [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb, T_ret)
#                #P_PV = curtailment*PV_out(T_amb, G*n_PV_modules)
#                P_BESS = BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)
#                [P_PV, P_BESS, P_Grid] = Electric_Power_Balance(P_Load, P_HP + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_BESS, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), 0, False, top)
#                [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1])                               
#                [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
                # print('No external control - Heating not required - Negative prices, charge BESS')
                # print('P_PV_av: ', P_PV_av)
    #            control_registry.append(5)
    
            elif (Cost_grid <= Cost_grid_Q1)*(not external_control):     # Low prices, charge BESS
                TESS_active = enable_TESS*False
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*(T_TESS_0 < T_TESS_max)*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4, Capacity_BESS_BoL, P_BESS_max, T_TESS_max)                
#                [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb, T_ret)
#                #P_PV = curtailment*PV_out(T_amb, G*n_PV_modules)
#                P_BESS = BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4
#                [P_PV, P_BESS, P_Grid] = Electric_Power_Balance(P_Load, P_HP + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_BESS, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), P_PV_av, False, top)
#                [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1])                               
#                [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
#                print('No external control - Heating not required - Low prices, charge BESS')
    #            control_registry.append(6)
    
            elif Cost_grid <= Cost_grid_Q3:     # Average prices, BESS compensate
                TESS_active = enable_TESS*False
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*False*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_BESS_ref, P_PV, setpoint] = check_P_BESS_setpoint(SoC_0, Capacity_BESS, P_Load, HP_Power(HP_active, T_amb, T_ret)[0], HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_PV_av, P_grid_DSO, P_BESS_max = P_BESS_max)
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, P_BESS_ref, Capacity_BESS_BoL, P_BESS_max, T_TESS_max, setpoint = True)        

#                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, 0, Capacity_BESS_BoL, P_BESS_max, T_TESS_max)
                # --- #
#                [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb, T_ret)
#                #P_PV = curtailment*PV_out(T_amb, G*n_PV_modules)
#                P_BESS = 0
#                [P_PV, P_BESS, P_Grid] = Electric_Power_Balance(P_Load, P_HP + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_BESS, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), P_PV_av, P_grid_DSO, top)
#                [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1])                               
#                [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)              
#                print('No external control - Heating not required - Average prices, BESS compensate')
                
    #            control_registry.append(7)
                
            else:                                   # High prices, discharge BESS
                TESS_active = enable_TESS*False
                HP_active = enable_HP*False
                HP_charge_TESS = enable_TESS*enable_HP*False*enable_HP_2_TESS
                [TESS_active, HP_active, HP_charge_TESS] = keep_minTin(T_in, TESS_active, HP_active, HP_charge_TESS, T_in_min = T_in_min)
                
                [P_HP, Qdot_HP, P_PV, P_BESS, P_Grid, T_TESS_new, Qdot_TESS, Qdot_TESS_SD, SoC_BESS, Capacity_BESS] = house_power_flow(enable_TESS, enable_HP, TESS_active, HP_active, HP_charge_TESS, T_amb, SoC_0, P_Load, P_PV_av, T_TESS_0, T_ret, T_soil, Capacity_BESS, P_grid_DSO, top, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4, Capacity_BESS_BoL, P_BESS_max, T_TESS_max)
#                [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb, T_ret)
#                #P_PV = curtailment*PV_out(T_amb, G*n_PV_modules)
#                P_BESS = BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max)/4
#                [P_PV, P_BESS, P_Grid] = Electric_Power_Balance(P_Load, P_HP + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_BESS, BESS_perm_max(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), BESS_perm_min(SoC_0, Capacity_BESS, P_BESS_max = P_BESS_max), P_PV_av, False, top)
#                [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1])                               
#                [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
#                print('No external control - Heating not required - High prices, discharge BESS')
#            control_registry.append(8)

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

        # [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb, T_ret)
        # [P_HP_TESS, Qdot_HP_TESS] = HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)
        
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
#    return [P_PV, P_PV_av, P_BESS, P_HP, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_Grid, SoC_BESS, Capacity_BESS, Qdot_HP, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1], T_TESS_new, Qdot_TESS, Qdot_TESS_SD, HP_active, HP_charge_TESS, TESS_active, Cost_grid, Cost_grid_median, Cost_grid_Q1, Cost_grid_Q3]



###############################################################################

def GA_control(T_in, T_set, T_amb, T_TESS_0, T_ret, T_soil, G, P_Load, P_grid_DSO, top, SoC_0, costs, Capacity_BESS, P_PV_av, horizon = 0, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP_2_TESS = True, external_control = False, SoC_BESS_min = 0.9, SoC_BESS_max = 0.2, m = 4000, TESS_min_op_T = 55 + 273, TESS_max_op_T = 90 + 273, T_TESS_min = 50 + 273, T_TESS_max = 90 + 273, m_dot = 0.22, c_f = 4200, dt = 0.25, House_type = 'apartment'):
    from MCES_library import HP_Power, update_TESS, new_house_Temperature, House_Thermal_Losses
    from GA_Optimization_lib import GA_Optimization
    
#    print('------GA control--------')
#    print([T_in, T_set, T_amb, T_TESS_0, T_ret, T_soil, G, P_Load, P_grid_DSO, top, SoC_0, costs, Capacity_BESS, P_PV_av])
    
#    print(SoC_0)
    
    best_candidate = GA_Optimization(T_amb, T_in, T_set, T_soil, P_PV_av, P_Load, G, T_TESS_0, SoC_0, horizon, costs) # consecutive_generations = 5, individuals = 200, beta = 6, theta_E = 0.483*0.25/2, theta_T = 1/1, theta_CO2 = 1/2.5)
    
    
    [HP_charge_TESS, TESS_active, HP_active, P_PV, P_BESS] = [best_candidate[0], best_candidate[1], best_candidate[2], best_candidate[3], best_candidate[4]]
    [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
    [P_HP, Qdot_HP] = HP_Power(HP_active, T_amb[0], T_ret)
    [T_TESS_new, Qdot_TESS, Qdot_TESS_SD] = update_TESS(TESS_active, T_ret, T_TESS_0, T_soil, 0, HP_Power(HP_charge_TESS, T_amb[0], T_TESS_0, T_sup = T_TESS_max)[1])
    T_in = new_house_Temperature(T_in, House_Thermal_Losses(T_in, T_amb[0], House_type = House_type )[0], House_Thermal_Losses(T_in, T_amb[0], House_type = House_type )[1], Qdot_HP = Qdot_HP, Qdot_TESS = Qdot_TESS)
    
#    print(P_BESS)
    
#    for i in [P_PV, P_BESS, P_HP, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0], P_Load + P_HP + HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[0] + P_BESS - P_PV, SoC_BESS, Capacity_BESS, House_Thermal_Losses(T_in, T_amb)[0], Qdot_HP, HP_Power(HP_charge_TESS, T_amb, T_TESS_0, T_sup = T_TESS_max)[1], T_TESS_new, Qdot_TESS, Qdot_TESS_SD, enable_HP_2_TESS, T_in]:
#        print(i)
    
    return [P_PV, P_BESS, P_HP, HP_Power(HP_charge_TESS, T_amb[0], T_TESS_0, T_sup = T_TESS_max)[0], P_Load[0] + P_HP + HP_Power(HP_charge_TESS, T_amb[0], T_TESS_0, T_sup = T_TESS_max)[0] + P_BESS - P_PV, SoC_BESS, Capacity_BESS, House_Thermal_Losses(T_in, T_amb[0], House_type = House_type)[0], Qdot_HP, HP_Power(HP_charge_TESS, T_amb[0], T_TESS_0, T_sup = T_TESS_max)[1], T_TESS_new, Qdot_TESS, Qdot_TESS_SD, enable_HP_2_TESS, T_in]
    


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
        
        return [-P_BESS, SoC_BESS, Capacity_BESS, P_Grid_sum, P_PV_av] # clip(P_PV_av, 0, -P_BESS)
        
#    [SoC_BESS, Capacity_BESS] = update_BESS_aeging(SoC_0, P_BESS, Capacity_BESS, Capacity_BESS_BoL)
        
###############################################################################
    
def Network_Control_Centralized(A, Z, I_0, V_0, selection, remaining, df_P_L, df_P_PV_av, df_P_PV, df_P_BESS, df_SOC_BESS, df_C_BESS, df_P_ref, timestep, t_0 = 0, Capacity_BESS_BoL = 1000, V_ref_min = 0.95*400, V_ref_max = 1.05*400):
    from pandas import DataFrame
    
    [V_n, I] = Network_state(A, Z, df_P_L, timestep, I_0, V_0)

    if (min(V_n) > V_ref_min) and (max(V_n) > V_ref_min):   # No action required
        for node in selection:        
        
            # Charge BESS
            [df_P_BESS.iloc[timestep,node], df_SOC_BESS.iloc[timestep,node], df_C_BESS.iloc[timestep,node], df_P_ref.iloc[timestep,node], df_P_PV.iloc[timestep,node]] = Central_BESS_Control(df_SOC_BESS.iloc[timestep-1,node], df_C_BESS.iloc[timestep-1,node], df_P_L.loc[timestep].sum(), df_P_PV_av.iloc[timestep,node], operation = "Charge")
        
    else:
        for node in selection:                
            #Discharge BESS
            [df_P_BESS.iloc[timestep,node], df_SOC_BESS.iloc[timestep,node], df_C_BESS.iloc[timestep,node], df_P_ref.iloc[timestep,node], df_P_PV.iloc[timestep,node]] = Central_BESS_Control(df_SOC_BESS.iloc[timestep-1,node], df_C_BESS.iloc[timestep-1,node], df_P_L.loc[timestep].sum(), df_P_PV_av.iloc[timestep,node], operation = "Support")

    
    return [df_P_BESS, df_SOC_BESS, df_C_BESS, df_P_ref, DataFrame(df_P_L.values + df_P_BESS.values, columns = df_P_L.columns)]

###############################################################################
#    
#def Network_Control(selection, df_P_Grid, df_P_n, df_V_n, df_I_n, df_P_L, df_P_ref, df_DSO_operation, ts, V_ref_min = 0.95, V_ref_max = 1.05, P_ref_min = -1.48, P_ref_max = 1.38):
def Network_Control(A, Z, I_0, V_0, Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control, P_BESS_max, Capacity_BESS_BoL, enable_HP, enable_TESS, V_ref_min = 0.95*400, V_ref_max = 1.05*400, P_ref_min = -1.48, P_ref_max = 1.38, force_temperature = False):
    from numpy import zeros
    B = zeros(len(I_0))
    B[0] = 1  
        
    [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in] = Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control = False, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, enable_HP = enable_HP, enable_TESS = enable_TESS, force_temperature = force_temperature)
    # Here is False for case 6
#    if False:
    if not follow_control:
#        print('Not controlled')
        return [[df_P_ref, df_DSO_operation], [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in]]
    
#    df_DSO_operation = assign_DF_values(df_DSO_operation, ts, selection, True)
#    df_P_ref = assign_DF_values(df_P_ref, ts, selection, 6)
    else: 
#        print('Controlled')
        [V_n, I] = Network_state(A, Z, df_P_Grid, timestep, I_0, V_0)
#        print('V_min = ', min(V_n))
#        print('V_max = ', max(V_n))
#        for node in range(len(df_P_Grid.loc[ts])):
#            if df_P_Grid.iloc[ts-1,node] > 0:
#                df_DSO_operation.iloc[ts,node] = True
#                df_P_ref.iloc[ts,node] = 6
#            else:
#                df_DSO_operation.iloc[ts,node] = False
#                df_P_ref.iloc[ts,node] = -5
        
        if (min(V_n) > V_ref_min) and (max(V_n) > V_ref_min):   # No action required
            # print('1')
            return [[df_P_ref, df_DSO_operation], [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in]]

        else:
            P_ref = Find_power_setpoints(A, Z, B, V_n, I, df_P_Grid.loc[timestep].tolist(), selection+remaining, selection)
            for node in selection:
                df_DSO_operation.iloc[timestep,node] = True
                df_P_ref.iloc[timestep,node] = P_ref[node-1]
            return [[df_P_ref, df_DSO_operation], Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, enable_HP = enable_HP, enable_TESS = enable_TESS, force_temperature = force_temperature)]        
#        elif (min(V_n) < V_ref_min) and (max(V_n) > V_ref_min):
#            # print('2')
#            P_ref = Find_power_setpoints(A, Z, B, V_n, I, df_P_Grid.loc[timestep].tolist(), selection+remaining, selection)
#            for node in selection:
#                if df_P_Grid.iloc[timestep,node] > P_ref_max:
#                    df_DSO_operation.iloc[timestep,node] = True
##                    df_P_ref.iloc[timestep,node] = P_ref_max   
#                    df_P_ref.iloc[timestep,node] = P_ref[node-1]
#            
#
#            return [[df_P_ref, df_DSO_operation], Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, enable_HP = enable_HP, enable_TESS = enable_TESS)]
#                    
#    #        return [df_P_ref, df_DSO_operation]    
#    
#        elif (min(V_n) > V_ref_min) and (max(V_n) < V_ref_min):
#            # print('3')
#            for node in selection:
#                if df_P_Grid.iloc[timestep,node] < P_ref_min:
#                    df_DSO_operation.iloc[timestep,node] = True
##                    df_P_ref.iloc[timestep,node] = P_ref_min
#                    df_P_ref.iloc[timestep,node] = P_ref[node-1]                    
#            return [[df_P_ref, df_DSO_operation], Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, enable_HP = enable_HP, enable_TESS = enable_TESS)]
#    #        return [df_P_ref, df_DSO_operation]    
#            
#        else:       #elif (df_V_n.loc[ts].min() > V_ref_min) and (df_V_n.loc[ts].max() > V_ref_min):    
#            # print('4')
#            for node in selection:
#                if df_P_Grid.iloc[timestep,node] > P_ref_max:
#                    df_DSO_operation.iloc[timestep,node] = True
##                    df_P_ref.iloc[timestep,node] = P_ref_max  
#                    df_P_ref.iloc[timestep,node] = P_ref[node-1]   
#                elif df_P_Grid.iloc[timestep,node] < P_ref_min:
#                    df_DSO_operation.iloc[timestep,node] = True
##                    df_P_ref.iloc[timestep,node] = P_ref_min
#                    df_P_ref.iloc[timestep,node] = P_ref[node-1]                
#            return [[df_P_ref, df_DSO_operation], Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, enable_HP = enable_HP, enable_TESS = enable_TESS)]
    #        return [df_P_ref, df_DSO_operation]                  
    #    return [df_P_ref, df_DSO_operation]

###############################################################################
def Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, ts, control_type = 'heuristic', follow_control = False, horizon = 0, energy_costs = [0, 0, 0, 0, 0.13, 0.483], CO2_costs = [0, 0, 0, 0, 0, 0.325], T_ret = 38 + 273, method = False, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP = True, enable_TESS = True, force_temperature = False, glazing_type = 'double', Cavity_Filling = False, Cavity_type = 'air', Wall_Cover = False, LRoof = 0.2, LCavity_wall = 0.05):
#    from MCES_library import House_Thermal_Losses
    
    if control_type == 'heuristic':
        if follow_control:
#            print('--- follow control ----')
            for node in range(len(df_T_in.loc[0])):
                if node in selection:
           
                    [df_P_PV.iloc[ts, node], df_P_BESS.iloc[ts, node], df_P_HP.iloc[ts, node], df_P_HP_TESS.iloc[ts, node], df_P_Grid.iloc[ts, node], df_SOC_BESS.iloc[ts, node], df_C_BESS.iloc[ts, node], df_Qdot_D.iloc[ts, node], df_Qdot_HP.iloc[ts, node], df_Qdot_HP_TESS.iloc[ts, node], df_T_TESS.iloc[ts, node], df_Qdot_TESS.iloc[ts, node], df_Qdot_TESS_SD.iloc[ts, node], df_Qdot_Boiler.iloc[ts, node], df_enable_HP_2_TESS.iloc[ts, node], df_T_in.iloc[ts, node]] = heuristic_control(df_T_in.iloc[ts-1, node], T_set[ts], T_amb[ts], df_T_TESS.iloc[ts-1, node], 38+273, T_soil[ts], 0, df_P_L.iloc[ts, node], df_P_ref.iloc[ts, node], df_DSO_operation.iloc[ts, node], df_SOC_BESS.iloc[ts-1, node], df_Prices['Current'][ts], df_Prices['Median'][ts], df_Prices['Q25'][ts], df_Prices['Q75'][ts], df_C_BESS.iloc[ts-1, node], df_P_PV_av.iloc[ts, node], enable_HP = enable_HP, enable_TESS = enable_TESS, enable_HP_2_TESS = df_enable_HP_2_TESS.iloc[ts-1, node], external_control = follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, House_type = df_House_type[Network_headers[5]][node+1], force_temperature = force_temperature)
    #               [P_PV,                      P_BESS,                 P_HP,               HP_Power(HP_charge_TESS)[0],    P_Grid,                     SoC_BESS,                   Capacity_BESS,      House_Thermal_Losses(T_in)[0], Qdot_HP,                     HP_Power(HP_charge_TESS)[1], T_TESS_new,               Qdot_TESS,                   Qdot_TESS_SD,                       Qdot_Boiler,            enable_HP_2_TESS,                       T_in]
                elif node in remaining:
                    [df_P_PV.iloc[ts, node], df_P_BESS.iloc[ts, node], df_P_HP.iloc[ts, node], df_P_HP_TESS.iloc[ts, node], df_P_Grid.iloc[ts, node], df_SOC_BESS.iloc[ts, node], df_C_BESS.iloc[ts, node], df_Qdot_D.iloc[ts, node], df_Qdot_HP.iloc[ts, node], df_Qdot_HP_TESS.iloc[ts, node], df_T_TESS.iloc[ts, node], df_Qdot_TESS.iloc[ts, node], df_Qdot_TESS_SD.iloc[ts, node], df_Qdot_Boiler.iloc[ts, node], df_enable_HP_2_TESS.iloc[ts, node], df_T_in.iloc[ts, node]] = [0, 0, 0, 0, df_P_L.iloc[ts, node], 0, 0, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1])[2], 0, 0, 0, 0, 0, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[0], False, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[1]]
    #                [Qdot_boiler, T_in_new, Qdot_D] = Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts])
                    
        else:
            for node in range(len(df_T_in.loc[0])):
                if node in selection:
           
                    [df_P_PV.iloc[ts, node], df_P_BESS.iloc[ts, node], df_P_HP.iloc[ts, node], df_P_HP_TESS.iloc[ts, node], df_P_Grid.iloc[ts, node], df_SOC_BESS.iloc[ts, node], df_C_BESS.iloc[ts, node], df_Qdot_D.iloc[ts, node], df_Qdot_HP.iloc[ts, node], df_Qdot_HP_TESS.iloc[ts, node], df_T_TESS.iloc[ts, node], df_Qdot_TESS.iloc[ts, node], df_Qdot_TESS_SD.iloc[ts, node], df_Qdot_Boiler.iloc[ts, node], df_enable_HP_2_TESS.iloc[ts, node], df_T_in.iloc[ts, node]] = heuristic_control(df_T_in.iloc[ts-1, node], T_set[ts], T_amb[ts], df_T_TESS.iloc[ts-1, node], 38+273, T_soil[ts], 0, df_P_L.iloc[ts, node], df_P_ref.iloc[ts, node], df_DSO_operation.iloc[ts, node], df_SOC_BESS.iloc[ts-1, node], df_Prices['Current'][ts], df_Prices['Median'][ts], df_Prices['Q25'][ts], df_Prices['Q75'][ts], df_C_BESS.iloc[ts-1, node], df_P_PV_av.iloc[ts, node], enable_HP = enable_HP, enable_TESS = enable_TESS, enable_HP_2_TESS = df_enable_HP_2_TESS.iloc[ts-1, node], external_control = follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, House_type = df_House_type[Network_headers[5]][node+1], force_temperature = force_temperature)
    #               [P_PV,                      P_BESS,                 P_HP,               HP_Power(HP_charge_TESS)[0],    P_Grid,                     SoC_BESS,                   Capacity_BESS,      House_Thermal_Losses(T_in)[0], Qdot_HP,                     HP_Power(HP_charge_TESS)[1], T_TESS_new,               Qdot_TESS,                   Qdot_TESS_SD,                       Qdot_Boiler,            enable_HP_2_TESS,                       T_in]
    
                elif node in remaining:
                    [df_P_PV.iloc[ts, node], df_P_BESS.iloc[ts, node], df_P_HP.iloc[ts, node], df_P_HP_TESS.iloc[ts, node], df_P_Grid.iloc[ts, node], df_SOC_BESS.iloc[ts, node], df_C_BESS.iloc[ts, node], df_Qdot_D.iloc[ts, node], df_Qdot_HP.iloc[ts, node], df_Qdot_HP_TESS.iloc[ts, node], df_T_TESS.iloc[ts, node], df_Qdot_TESS.iloc[ts, node], df_Qdot_TESS_SD.iloc[ts, node], df_Qdot_Boiler.iloc[ts, node], df_enable_HP_2_TESS.iloc[ts, node], df_T_in.iloc[ts, node]] = [0, 0, 0, 0, df_P_L.iloc[ts, node], 0, 0, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1])[2], 0, 0, 0, 0, 0, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[0], False, Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts], House_type = df_House_type[Network_headers[5]][node+1], glazing_type = glazing_type, Cavity_Filling = Cavity_Filling, Cavity_type = Cavity_type, Wall_Cover = Wall_Cover, LRoof = LRoof, LCavity_wall = LCavity_wall)[1]]
                    
    #                [Qdot_boiler, T_in_new, Qdot_D] = Gas_Boiler(df_T_in.iloc[ts-1, node], T_amb[ts], T_set[ts])
                
#        return [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in]


    elif control_type == 'GA':
        from GA_Optimization_lib import TS_forecast
        
        h = ts + horizon + 1        
        
        horizon_costs = [[],[]]
        for i in range(horizon + 1):
            energy_costs[-1] = df_Prices['Current'][ts + horizon]
            horizon_costs[0].append(energy_costs)        
            horizon_costs[1].append(CO2_costs)
            
            
            for node in range(len(df_T_in.loc[0])):            
                if node in selection:    
                    
#                    print("____node_____")
#                    print(df_T_in.iloc[ts-1, node])
#                    print(T_set[ts])
#                    print(T_amb[ts])
#                    print(df_T_TESS.iloc[ts-1, node])
#                    print(T_ret)
#                    print(T_soil[ts])
#                    print(df_P_L.iloc[ts, node])
#                    print(df_P_ref.iloc[ts, node])
#                    print(df_DSO_operation.iloc[ts, node])
#                    print(df_SOC_BESS.iloc[ts-1, node])
#                    print(horizon_costs)
#                    print(df_C_BESS.iloc[ts-1, node])
#                    print(df_P_PV_av.iloc[ts, node])
#                    print(TS_forecast(T_amb[ts:h], ts, method, forecast_type = 'Temperature'))
#                    print(df_P_PV_av.iloc[ts:h, node].tolist())
#                    print(G[ts:h])
#                    print(df_P_L.iloc[ts:h, node].tolist())
#                    
#                    print("____node_____")                    
#                    GA_control(df_T_in.iloc[ts-1, node], T_set[ts], T_amb[ts], df_T_TESS.iloc[ts-1, node], T_ret, T_soil[ts], 0, df_P_L.iloc[ts, node], df_P_ref.iloc[ts, node], df_DSO_operation.iloc[ts, node], df_SOC_BESS.iloc[ts-1, node], horizon_costs, df_C_BESS.iloc[ts-1, node], df_P_PV_av.iloc[ts, node])
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
#                   [P_PV,                      P_BESS,                 P_HP,                   HP_Power(HP_charge_TESS)[0], P_Load +,                  SoC_BESS,                   Capacity_BESS,          House_Thermal_Losses()[0], Qdot_HP,                     HP_Power(HP_charge_TESS)[1], T_TESS_new,                Qdot_TESS,                  Qdot_TESS_SD,                   enable_HP_2_TESS,                       T_in]              = GA_control(T_in,                       T_set,      T_amb,  T_TESS_0,                     T_ret, T_soil,    G, P_Load,                  P_grid_DSO,         top,                        SoC_0,                      costs,          Capacity_BESS,              P_PV_av, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP_2_TESS = True, external_control = False, SoC_BESS_min = 0.9, SoC_BESS_max = 0.2, m = 4000, TESS_min_op_T = 55 + 273, TESS_max_op_T = 90 + 273, T_TESS_min = 50 + 273, T_TESS_max = 90 + 273, m_dot = 0.22, c_f = 4200, dt = 0.25):

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
    
#                                                        [0,         0,          0,          0,         0,  False,              0.5,        10,         0,      0,              True,               75+273, 20+273+
    
    df_P_Grid.loc[0] = df_P_L.loc[0]
    df_P_ref = assign_DF_values(df_P_ref, 0, selection, 1.5)#df_P_L.loc[0]
    
    for node in selection:
#        House_state = House_Thermal_Losses(df_T_in.iloc[0, node], T_amb[0])
        df_Qdot_D.iloc[0, node] = House_Thermal_Losses(df_T_in.iloc[0, node], T_amb[0], House_type = df_House_type[Network_headers[5]][node+1])[0]#House_state[0]
        df_Qdot_TESS_SD.iloc[0, node] = Qdot_SD_TESS2Soil(T_soil[0], initial_conditions[10] - 273)

    for node in remaining:
        df_T_in.iloc[0, node] = 19+273
        df_Qdot_D.iloc[0, node] = House_Thermal_Losses(df_T_in.iloc[0, node], T_amb[0], House_type = df_House_type[Network_headers[5]][node+1])[0]#House_state[0]
    
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
def simulate_Network(DF_Network, Network_headers, Wires_headers, selection, remaining, S_n, t_simulation, V_0 = False, I_0 = False, t_0 = 0, wire_file = 'Wires.xlsx', H = 0, control_type = 'heuristic', follow_control = False, P_BESS_max = 10, Capacity_BESS_BoL = 10, enable_HP = True, enable_TESS = True, V_feeder = 400, dt = 0.25, force_temperature = False):
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
#    df_P_PV_av = Create_PV_DF(S_n)
    df_P_PV = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)    
    df_P_L = DataFrame(S_n[t_0:t_0+t_simulation+H], columns = labels)#.div(1000)
#    df_P_L = DataFrame(S_n, columns = labels)#.div(1000)
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
    
#    df_House_type = DF_Network[[Network_headers[1], Network_headers[5]]].sort_values(by=[Network_headers[1]])
#    df_House_type = df_House_type[df_House_type[Network_headers[5]].notnull()].set_index(Network_headers[1])
    
    df_House_type = DF_Network[[Network_headers[1], Network_headers[5]]].sort_values(by=[Network_headers[1]])
    df_House_type = df_House_type.set_index(Network_headers[1])
    
    [df_P_HP, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_T_in] = set_initialConditions(Network_headers, selection, remaining, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_T_in, T_amb, T_soil, df_House_type)
    
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
    
    P_ref_min = -1.5
    P_ref_max = 1.5
        
    print('Start Simulation')   
#    if control_type == 'GA':
#        global population_registry
#        population_registry = []     
#        population_registry.append([])            

    for timestep in range(1,int(t_simulation)):
#        print('---', timestep, '---') # , 4*df_Prices['Current'][timestep]
        if (timestep)%(24/dt) == 0:
            print('Day', t_0/(24/dt) + int(timestep/(24/dt)))
        start_ts = time.time()
        
#        [df_P_ref, df_DSO_operation] = Network_Control(selection, df_P_Grid, df_P_n, df_V_n, df_I_n, df_P_L, df_P_ref, df_DSO_operation, timestep)
        
        
#        [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in] = Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, enable_HP = enable_HP, enable_TESS = enable_TESS)
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

def simulate_Network_Centralized_BESS(DF_Network, Network_headers, Wires_headers, selection, remaining, S_n, P_PV_peak, df_P_L, t_simulation, V_0 = False, I_0 = False, t_0 = 0, wire_file = 'Wires.xlsx', H = 0, control_type = 'heuristic', follow_control = False, P_BESS_max = 1000, Capacity_BESS_BoL = 1000, V_feeder = 400, dt = 0.25):
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
        

        [df_P_BESS, df_SOC_BESS, df_C_BESS, df_P_ref, df_P_Grid] = Network_Control_Centralized(A, Z, I_0, V_0, selection, remaining, df_P_L, df_P_PV_av, df_P_PV, df_P_BESS, df_SOC_BESS, df_C_BESS, df_P_ref, timestep, t_0, Capacity_BESS_BoL = 1000, V_ref_min = 0.95*400, V_ref_max = 1.05*400)
        
        
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

def Simulate_Network_Sizing(DF_Network, Network_headers, Wires_headers, S_n, t_simulation, V_0 = False, I_0 = False, t_0 = 0, wire_file = 'Gaia_cables.xlsx', P_BESS_max = 10, Capacity_BESS_BoL = 10, V_feeder = 400, dt = 0.25):

    from numpy import empty, zeros
    from pandas import DataFrame
    import csvreader
    
#    print('Creating network')    
    
    [A, Z] = Create_network(DF_Network, Network_headers, Wires_headers, wire_file)


#    print('Loading data')
    CSVDataPrices = csvreader.read_data(csv='DA_Prices_15min.csv', address='')
    CSVDataPrices.data2array()
    Energy_price = [i[0]*dt for i in CSVDataPrices.ar[t_0:t_0+t_simulation]] # [energy_costs[-1]*0.25 for i in CSVDataPrices.ar]

#    print('Creating dataframes')   
    
    labels = [str(i+1) for i in range(len(S_n[0]))]

    # Only for nodes with load
    df_P_PV_av = Create_PV_DF(S_n, t_0, t_simulation)
    df_P_PV = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)    
#    df_P_L = DataFrame(S_n[t_0:t_0+t_simulation+H], columns = labels)#.div(1000)
    df_P_BESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
#    df_P_Grid = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_P_Grid = DataFrame(S_n[t_0:t_0+t_simulation], columns = labels)#.div(1000)    
    df_SOC_BESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_C_BESS = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)

    
    # For every node
    df_P_n = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    df_V_n = DataFrame(zeros([int(t_simulation), len(A)]))
    df_I_n = DataFrame(zeros([int(t_simulation), len(A)]), columns = labels)
    
    df_Prices = create_pricesDF(Energy_price, t_simulation)    

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
         
        
    print('Start Simulation')          

    for timestep in range(1,int(t_simulation)):
#        print('---', timestep, '---')
        if (timestep)%(24/dt) == 0:
            print('Day', t_0/(24/dt) + int(timestep/(24/dt)))
        start_ts = time.time()
        
#        [df_P_ref, df_DSO_operation] = Network_Control(selection, df_P_Grid, df_P_n, df_V_n, df_I_n, df_P_L, df_P_ref, df_DSO_operation, timestep)                
#        [df_P_PV, df_P_BESS, df_P_HP, df_P_HP_TESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_T_TESS, df_Qdot_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in] = Local_Control(Network_headers, selection, remaining, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_enable_HP_2_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_T_in, df_Prices, T_amb, G, T_set, T_soil, df_House_type, timestep, control_type, follow_control, P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, enable_HP = enable_HP, enable_TESS = enable_TESS)
        
        [V_n, I] = Network_state(A, Z, df_P_Grid, timestep, I_0, V_0)
        
        V_n_registry[:,timestep] = V_n
        I_registry[:,timestep] = I
        
        V_0 = V_n
        I_0 = I
        
        end_ts = time.time()
        t_registry.append(end_ts-start_ts)

    df_V_n = DataFrame(V_n_registry).T
    df_V_n.columns = labels
    df_I_n = DataFrame(I_registry).T
    df_I_n.columns = labels
    
    return [V_n_registry, I_registry, df_P_PV_av, df_P_PV, df_P_BESS, df_P_Grid, df_SOC_BESS, df_C_BESS, df_P_n, df_V_n, df_I_n, df_Prices]

###############################################################################

def Voltage_heatmap(V_n_registry, V_feeder = 400):
    
    from numpy import array
    import matplotlib.pyplot as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })      
    
    V_n_registry = array(V_n_registry)/V_feeder - 1
    
    plt.figure(constrained_layout=True)
    plt.imshow(V_n_registry)
    plt.colorbar()
    plt.show()

###############################################################################
def Voltage_Power_models(Case_S, Case_W, Case_Y, all_nodes, V_feeder = 400):
    
    from numpy import linspace, poly1d, polyfit
    import matplotlib.pyplot as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })      
    
#    P_av_C0 = linspace(0, 3, 300)
#    P_av_C1 = linspace(-6, 3, 900)
#    P_av_C2 = linspace(-6, 8, 1400)
#    P_av_C3 = linspace(-6, 8, 1400)
#    P_av_C5 = linspace(-4.5, 8.5, 1300)
#    
##    P_av = linspace(-6, 8.5, 1450)
#    
#    V_C0 = [-0.0007060*P**2 - 0.03396*P + 0.9997 for P in P_av_C0]
#    V_C1 = [-0.0007895*P**2 - 0.03373*P + 0.9997 for P in P_av_C1]
#    V_C2 = [-0.0007996*P**2 - 0.03373*P + 0.9996 for P in P_av_C2]
#    V_C3 = [-0.0008002*P**2 - 0.03371*P + 0.9996 for P in P_av_C3]
#    V_C5 = [-0.0010600*P**2 - 0.03379*P + 1 for P in P_av_C5]
#    
#
#    plt.figure(constrained_layout=True)
#    plt.grid()
#    plt.plot(P_av_C0, V_C0, label='Base', color ='k')
#    plt.plot(P_av_C1, V_C1, label='Scenario 1', color ='tab:blue')
#    plt.plot(P_av_C2, V_C2, label='Scenario 2', color ='tab:orange')
#    plt.plot(P_av_C3, V_C3, label='Scenario 3', color ='tab:red')
#    plt.plot(P_av_C5, V_C5, label='Scenario 5', color ='tab:gray')
#    plt.xlabel('Average grid power, $\overline{P}_{Grid}$, [kW]')
#    plt.ylabel('Worst voltage, $\hat{V}$, [V]')
#    plt.legend(loc='lower left')
##    plt.title('Model: '+"${}$".format(eq_latex))
#    plt.show()    


#    P_av = linspace(-6, 8.5, 1450)
    
    plt.figure(constrained_layout=True)
    plt.grid()    

    
    for scenario in Case_Y:

        l_P = []
        l_V = []
        
        for P,V in zip(scenario[12][[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), [j/V_feeder for j in scenario[27][[str(i + 1) for i in all_nodes]].min(axis=1).tolist()]):
            if P > 0:        
                l_P.append(P)
                l_V.append(V)
        for P,V in zip(scenario[12][[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), [j/V_feeder for j in scenario[27][[str(i + 1) for i in all_nodes]].max(axis=1).tolist()]):
            if P < 0:        
                l_P.append(P)
                l_V.append(V)    
            
        model = poly1d(polyfit(l_P, l_V, 2))
        
#        print(model)        
        plt.plot(linspace(min(l_P), max(l_P), 130), model(linspace(min(l_P), max(l_P), 130)), color ='k', label='Year') # , label=eq_latex
        
    
    for scenario in Case_S:

        l_P = []
        l_V = []
        
        for P,V in zip(scenario[12][[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), [j/V_feeder for j in scenario[27][[str(i + 1) for i in all_nodes]].min(axis=1).tolist()]):
            if P > 0:        
                l_P.append(P)
                l_V.append(V)
        for P,V in zip(scenario[12][[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), [j/V_feeder for j in scenario[27][[str(i + 1) for i in all_nodes]].max(axis=1).tolist()]):
            if P < 0:        
                l_P.append(P)
                l_V.append(V)    
            
        model = poly1d(polyfit(l_P, l_V, 2))
        
#        print(model)        
        plt.plot(linspace(min(l_P), max(l_P), 130), model(linspace(min(l_P), max(l_P), 130)), color ='r', alpha = 0.5) # , label=eq_latex
        
    plt.plot(0,1, color ='r', label='Summer')

    for scenario in Case_W:

        l_P = []
        l_V = []
        
        for P,V in zip(scenario[12][[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), [j/V_feeder for j in scenario[27][[str(i + 1) for i in all_nodes]].min(axis=1).tolist()]):
            if P > 0:        
                l_P.append(P)
                l_V.append(V)
        for P,V in zip(scenario[12][[str(i + 1) for i in all_nodes]].mean(axis=1).tolist(), [j/V_feeder for j in scenario[27][[str(i + 1) for i in all_nodes]].max(axis=1).tolist()]):
            if P < 0:        
                l_P.append(P)
                l_V.append(V)    
            
        model = poly1d(polyfit(l_P, l_V, 2))
        
#        print(model)        
        plt.plot(linspace(min(l_P), max(l_P), 130), model(linspace(min(l_P), max(l_P), 130)), color ='b', alpha = 0.5) # , label=eq_latex

    plt.plot(0,1, color ='b', label='Winter')
        
    plt.xlabel('Average grid power, $\overline{P}_{Grid}$, [kW]')
    plt.ylabel('Worst voltage, $\hat{V}$, [V]')
    plt.legend(loc='lower left')
#    plt.title('Model: '+"${}$".format(eq_latex))
    plt.show()        
    
###############################################################################    

def Voltage_Power_correlation(all_nodes, df_P_Grid, df_V_n, df_Prices, V_feeder = 400, t0 = 0, tf = 1*24*4):
    
#    from sympy import S, symbols, printing
    from numpy import linspace, poly1d, polyfit
    import matplotlib.pyplot as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })      

    fig, axs = plt.subplots(3, 1, sharex=True) # , constrained_layout=True
    
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
    
    print(model)
    
    
#    p = polyfit(l_P, l_V, 2)
#    x = symbols("x")
#    poly = sum(S("{:6.2f}".format(v))*x**i for i, v in enumerate(p[::-1]))
#    eq_latex = printing.latex(poly)
    
    
        
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
    
    
    
#    axs[0, 0].plot(P_L) # Load
#    axs[0, 0].set_ylabel('$P_{L}$, [kW]')
#    axs[1, 0].plot(P_PV_av, color='r', label='av') # PV av
#    axs[1, 0].plot(P_PV, color='b', label='PV') # PV av, PV
#    axs[1, 0].set_ylabel('$P_{PV}$, [kW]')
#    axs[1, 0].legend(loc='upper right', ncol = 1, fontsize=10)
#    axs[2, 0].plot(P_HP, color='r', label='L') # HP
#    axs[2, 0].plot(P_HP_TESS, color='b', label='TESS') # TESS    
#    axs[2, 0].set_ylabel('$P_{HP}$, [kW]')
#    axs[2, 0].legend(loc='upper right', ncol = 1, fontsize=10)
#    axs[3, 0].plot(P_BESS) # BESS
#    axs[3, 0].set_ylabel('$P_{BESS}$, [kW]')
#    axs[4, 0].plot(SoC_BESS) # BESS
#    axs[4, 0].set_ylabel('$SoC_{BESS}$, [%]')
#    axs[4, 0].set_ylim([15,95])    
#    axs[5, 0].plot(P_Grid)
#    axs[5, 0].set_ylabel('$P_{Grid}$, [kW]')
#    axs[6, 0].plot(Prices)
#    axs[6, 0].set_ylabel('c, [€/kWh]')
#    axs[6, 0].set_xlim([start,end])    
    
###############################################################################

#Voltage_Power_models(Case_4_MCES_NA_S, Case_4_MCES_NA_W, Case_4_MCES_NA_Y, all_nodes, V_feeder = 400)
#Voltage_Power_correlation(all_nodes, Case_4_MCES_NA_Y[0][12], Case_4_MCES_NA_Y[0][27], Case_4_MCES_NA_Y[0][29], tf = t_final)

###############################################################################

def Simulate_loads(V_0, S_ref):
    
    # Include here the code to estimate the outputs of each building.
    # A recursive function could be used, or even the simulation of the
    # non-controlled network for the control.
    
    return False

###############################################################################

def Network_control(V_0):
    
    # Include here the code to estimate the desired power per building, based
    # on the current status of the network, and maybe other parameters.
    return False

###############################################################################
    
def simulate_Controlled_Network(A, Z, t_simulation, V_0 = False, I_0 = False, t_0 = 0, dt = 0.25):

    from numpy import empty
    
    V_n_registry = empty([len(A), int(t_simulation)])
    I_registry = empty([len(A), int(t_simulation)])      
    
    for timestep in range(int(t_simulation)):

        S_ref = Network_control(V_0)        
        
        S_n = Simulate_loads(V_0, S_ref)
        
        [V_n, I] = Network_state(A, Z, S_n, timestep + t_0, I_0, V_0)
        
        V_n_registry[:,timestep] = V_n
        I_registry[:,timestep] = I    
    
        return [V_n_registry, I_registry]

###############################################################################
        
def Create_Load_Profiles(file_name):
    import csvreader  
    PowerData = csvreader.read_data(csv=file_name, address='')
    PowerData.data2array()   
    
    return PowerData.ar*1000

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

def Save_CSV(variable, file_name = ''):
    import csv

    with open(file_name, 'w') as file:
        writer = csv.writer(file)
        writer.writerows(variable)    

###############################################################################

def create_random_selection(S_n, percentage = 100):
    from random import sample
    from math import ceil
    
    Total_nodes = []
    
    for node in range(len(S_n[0])):
        if max(S_n[:,node]) != 0:
            Total_nodes.append(node)
    
    
    return sorted(sample(Total_nodes, ceil(percentage*len(Total_nodes)/100)))

###############################################################################

def create_device_list(S_n, selection, value = 1):
    from numpy import zeros
    
    device_list = zeros(len(S_n[0]))
    for node in selection:
        device_list[node] = value
    
    return device_list

###############################################################################

def create_base_selections(S_n, p_ranges = 10, n_samples = 5, csv_create = False, csv_labels = ['Stedin_301_all_nodes.csv', 'Stedin_301_selected_nodes_list.csv', 'Stedin_301_remaining_nodes_list.csv']):
    
    penetrations = [10*i for i in range(1, p_ranges)]

        
    all_nodes = create_random_selection(S_n, 100)
    selected_nodes_list = [[create_random_selection(S_n, penetration) for sample in range(n_samples)] for  penetration in penetrations]
    remaining_nodes_list = [[[node for node in all_nodes if node not in sample] for sample in penetration] for penetration in selected_nodes_list]
    
    if csv_create:
        Save_CSV([all_nodes], file_name = csv_labels[0])
        Save_CSV([selected_nodes_list], file_name = csv_labels[1])
        Save_CSV([remaining_nodes_list], file_name = csv_labels[2])
        
    
    return [all_nodes, selected_nodes_list, remaining_nodes_list]

###############################################################################

def select_nodes(penetration, case, selected_nodes_list, remaining_nodes_list, p_ranges = 10):
    
    penetrations = [10*i for i in range(1, p_ranges)]
    
    return [selected_nodes_list[penetrations.index(penetration)][case], remaining_nodes_list[penetrations.index(penetration)][case]]
    

###############################################################################
def load_base_selections():
    import csvreader
    all_nodes_Data = csvreader.read_data(csv = 'Stedin_301_all_nodes.csv', address = '', delim = ',')
    all_nodes_Data.data2array()
    
    selected_nodes_Data = csvreader.read_data(csv = 'Stedin_301_selected_nodes_list.csv', address = '', delim = ',')
    selected_nodes_Data.data2array()
    
    remaining_nodes_Data = csvreader.read_data(csv = 'Stedin_301_remaining_nodes_list.csv', address = '', delim = ',')
    remaining_nodes_Data.data2array()    
    
    
#    all_nodes_list = all_nodes_Data.ar.tolist()[0]
#    selected_nodes_list = selected_nodes_Data.ar.tolist()[0]
#    remaining_nodes_list = remaining_nodes_Data.ar.tolist()[0]
    
    return [all_nodes_Data.ar.tolist()[0], selected_nodes_Data.ar.tolist()[0], remaining_nodes_Data.ar.tolist()[0]]
    
###############################################################################

def create_histograms_consumption(csv_name = 'consumptions.csv', csv_address = '', csv_delim = ',', base=10, return_val = True):
    import csvreader    
    from numpy import arange, linspace
    import matplotlib.pylab as plt
    
    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    }) 
    
    Consumption_Data = csvreader.read_data(csv = csv_name, address = csv_address, delim = csv_delim)
    Consumption_Data.data2cols()
    Consumptions = [i/1000 for i in Consumption_Data.col[0]]
    
    bins = linspace(0, int(max(Consumptions))+1, 33)
    
    plt.figure(constrained_layout=True)    
    plt.hist(Consumptions, bins, alpha = 0.5, density = False, histtype = "bar", cumulative = False)
    plt.xticks(arange(0, max(Consumptions)+2, step = 1))
    plt.xlim((0, max(Consumptions)+1))
    plt.yticks(arange(0, 18+3, step = 2))
    plt.ylim((0, 20))
    plt.xlabel('Yearly consumption [MWh]')
    plt.ylabel('Frequency')
    plt.show()    
    
    if return_val:
        return Consumptions
    
###############################################################################

def create_histograms_samples(penetration, all_nodes, selected_nodes_list, p_ranges = 10, base=10):

    from numpy import arange, linspace
    import matplotlib.pylab as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })      
    
    penetrations = [10*i for i in range(1, p_ranges)]        
    samples = [[all_nodes.index(node) for node in selected_nodes] for selected_nodes in selected_nodes_list[penetrations.index(penetration)]]

    bins = linspace(0, base*round(len(all_nodes)/base), round(len(all_nodes)/base))
    colors = ['black', 'coral', 'gold', 'silver', 'blue']    
    
    plt.figure(constrained_layout=True)
#    for sample in samples:        
    plt.hist(samples, bins, alpha=0.5, density=False, histtype="bar", cumulative=False, color=colors, label=['a', 'b', 'c', 'd', 'e'])
    plt.xticks(arange(0, base*round(len(all_nodes)/base), step=base))
    plt.xlim((0, base*round(len(all_nodes)/base)))
    plt.legend(['a', 'b', 'c', 'd', 'e'], loc='lower center', ncol=5)
    plt.title('Penetration percentage '+str(penetration))
    plt.xlabel('Node number')
    plt.ylabel('Samples')
    plt.show()
    
#    return l

###############################################################################    


def p_value_lists(case_data, cases = 5, penetrations = 4, seasons = 2, index=4, plot = False):
    from scipy.stats import f_oneway
    
    p_value_list=[]
    
    for season in range(seasons):
        for penetration in range(penetrations):
#            print(season*penetrations*cases + penetration*cases + 0, season*penetrations*cases + penetration*cases + 1, season*penetrations*cases + penetration*cases + 2, season*penetrations*cases + penetration*cases + 3, season*penetrations*cases + penetration*cases + 4 )
            F, p = f_oneway(case_data[season*penetrations*cases + penetration*cases + 0][index],case_data[season*penetrations*cases + penetration*cases + 1][index],case_data[season*penetrations*cases + penetration*cases + 2][index],case_data[season*penetrations*cases + penetration*cases + 3][index],case_data[season*penetrations*cases + penetration*cases + 4][index])
            p_value_list.append(100*len([i for i in p if i>=0.05])/len([i for i in F if i>=0]))

            if plot:            
                boxplot_headers = ['a', 'b','c', 'd', 'e']
                Make_voltage_boxplots([case_data[season*penetrations*cases + penetration*cases + 0][index],case_data[season*penetrations*cases + penetration*cases + 1][index],case_data[season*penetrations*cases + penetration*cases + 2][index],case_data[season*penetrations*cases + penetration*cases + 3][index],case_data[season*penetrations*cases + penetration*cases + 4][index]], boxplot_headers)
            
    return p_value_list

###############################################################################

def case_parameters(case):
    if case == 1 or case == 0:
        enable_HP = False
        BESS_power = 0      # 10
        BESS_energy = 10
        enable_TESS = False
        TESS_capacity = 4
        follow_control = False
    elif case == 2:
        enable_HP = True
        BESS_power = 0      # 10
        BESS_energy = 10
        enable_TESS = False
        TESS_capacity = 4     
        follow_control = False
    elif case == 3:
        enable_HP = True
        BESS_power = 10      # 10
        BESS_energy = 10
        enable_TESS = False
        TESS_capacity = 4 
        follow_control = False
    elif case == 4:
        enable_HP = True
        BESS_power = 10      # 10
        BESS_energy = 10
        enable_TESS = True
        TESS_capacity = 4
        follow_control = False
    elif case == 5:
        enable_HP = True
        BESS_power = 10      # 10
        BESS_energy = 10
        enable_TESS = True
        TESS_capacity = 4
        follow_control = True   
    elif case == 6:
        enable_HP = True
        BESS_power = 10      # 10
        BESS_energy = 10
        enable_TESS = True
        TESS_capacity = 4
        follow_control = False        
    elif case == 7:
        enable_HP = True
        BESS_power = 10      # 10
        BESS_energy = 10
        enable_TESS = False
        TESS_capacity = 4 
        follow_control = True    
        
        
    elif case == 'CIGRE_18':
        enable_HP = True
        BESS_power = 10      # 10
        BESS_energy = 10
        enable_TESS = True
        TESS_capacity = 4 
        follow_control = False        

    return [enable_HP, BESS_power, BESS_energy, enable_TESS, TESS_capacity, follow_control]  
     
###############################################################################
    
def compare_ref(case_n, case_0):
    
    return [[100*c/r for c,r in zip(penetration, ref)] for penetration, ref in zip(case_n, case_0)]

###############################################################################
def count_V_compliance(Case_n, lim = 0.05, percentage = True, index = 4, V_feeder = 400, verbose = False):
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

    if verbose:
        print('Voltage ', 1-lim)
        print('Penetration: 20% ', min(v_low[0:5]), max(v_low[0:5]), len(v_low[0:5]))
        print('Penetration: 40% ', min(v_low[5:10]), max(v_low[5:10]), len(v_low[5:10]))
        print('Penetration: 60% ', min(v_low[10:15]), max(v_low[10:15]), len(v_low[10:15]))
        print('Penetration: 80% ', min(v_low[15:20]), max(v_low[15:20]), len(v_low[15:20]))
    
    
        print('Voltage ', 1+lim)
        print('Penetration: 20% ', min(v_up[0:5]), max(v_up[0:5]), len(v_up[0:5]))
        print('Penetration: 40% ', min(v_up[5:10]), max(v_up[5:10]), len(v_up[5:10]))
        print('Penetration: 60% ', min(v_up[10:15]), max(v_up[10:15]), len(v_up[10:15]))
        print('Penetration: 80% ', min(v_up[15:20]), max(v_up[15:20]), len(v_up[15:20]))

    return [v_low, v_up]
###############################################################################

#Performance_comparison(Case_n, Case_0_base_W, Case_1_PV_W, Case_2_PV_HP_W, Case_3_PV_HP_BESS_W, Case_4_MCES_NA_W, Case_5_MCES_PA_W, Case_6_MCES_UA_W, Case_7_PV_HP_BESS_A_S)

def Performance_comparison(Case_n, Case_0, Case_1 = False, Case_2 = False, Case_3 = False, Case_4 = False, Case_5 = False, Case_6 = False, penetrations = 4, cases = 5, indexes = [4,4, 25, 12, 12, 12, 7, 11, 16, 9, 10, 19], dt = 0.25, verbose = True, plots = True, gas_price = 1.44):
    
    
    [V_95, V_105] = count_V_compliance(Case_n, 0.05, verbose = verbose)
    [V_90, V_110] = count_V_compliance(Case_n, 0.1, verbose = verbose)
    
    E_prices = Case_n[0][29]['Current'].tolist()
    E_load_ref = [(365/7)*dt*sum([i for i in Case_0[0][8][str(node+1)].tolist()]) for node in Case_n[-1][3]]            

    if not Case_1:
        gas_price = 0
        
    precharge_TESS = 0
#    if Case_3:
#        precharge_TESS = (4000*4186*25/3600000)*Case_n[0][29]['Current'].mean()/dt
#    else:
#        precharge_TESS = 0
    
    P_cable = []
    
    Delta_Cost = []
    Delta_Cost_gas = []
    Delta_E_grid_in = []
    Delta_E_grid_out = []
    Delta_E_PV_curtailed = []
    Delta_E_HP = []
    Delta_E_BESS_in = []
    Delta_C_BESS = []
    Delta_E_HP_TESS = []
    Delta_Qdot_TESS = []
    

    # note that, when the reference is not case 0, the first index of Case_0 is not [0] but [penetration*penetrations + penetration + case]    
    for penetration in range(penetrations):
        for case in range(cases):
#            E_grid_ref.append([dt*sum(Case_0[0][12][str(node+1)].tolist()) for node in Case_n[penetration*penetrations + penetration + case][3]])
            
#            P_grid.append([Case_n[penetration*penetrations + penetration + case][12][str(node+1)].tolist() for node in Case_n[penetration*penetrations + penetration + case][3]])
#            P_grid_ref.append([Case_0[0][12][str(node+1)].tolist() for node in Case_n[penetration*penetrations + penetration + case][3]])
            P_grid_ref = [Case_0[0][12][str(node+1)].tolist() for node in Case_n[penetration*penetrations + penetration + case][3]]
            Qdot_boiler_ref = [Case_0[0][23][str(node+1)].tolist() for node in Case_n[penetration*penetrations + penetration + case][3]]
            P_grid = [Case_n[penetration*penetrations + penetration + case][12][str(node+1)].tolist() for node in Case_n[penetration*penetrations + penetration + case][3]]            
            
            C_grid_ref = [sum([dt*p*c for p,c in zip(profile, E_prices)]) for profile in P_grid_ref]
            C_grid_ref_gas = [sum([dt*(p*c + gas_price*q/(1000*12.7)) for p,c,q in zip(profile, E_prices, profile_gas)]) for profile,profile_gas in zip(P_grid_ref, Qdot_boiler_ref)]
            C_grid = [sum([dt*p*c for p,c in zip(profile, E_prices)])+precharge_TESS for profile in P_grid]

#            E_grid.append([dt*sum(Case_n[penetration*penetrations + penetration + case][12][str(node+1)].tolist()) for node in Case_n[penetration*penetrations + penetration + case][3]])
#            E_grid_in.append([dt*sum([i for i in Case_n[penetration*penetrations + penetration + case][12][str(node+1)].tolist() if i>0]) for node in Case_n[penetration*penetrations + penetration + case][3]])
            
            E_grid_ref = [dt*sum(Case_0[0][12][str(node+1)].tolist()) for node in Case_n[penetration*penetrations + penetration + case][3]]
            E_grid = [dt*sum(Case_n[penetration*penetrations + penetration + case][12][str(node+1)].tolist()) for node in Case_n[penetration*penetrations + penetration + case][3]]
            
            E_grid_in_ref = [dt*sum([i for i in Case_0[0][12][str(node+1)].tolist() if i>0]) for node in Case_n[penetration*penetrations + penetration + case][3]]            
            E_grid_in = [dt*sum([i for i in Case_n[penetration*penetrations + penetration + case][12][str(node+1)].tolist() if i>0]) for node in Case_n[penetration*penetrations + penetration + case][3]]

            Delta_Cost.append([100*(c/c_ref-1) for c,c_ref in zip(C_grid, C_grid_ref)])          
            Delta_Cost_gas.append([100*(c/c_ref-1) for c,c_ref in zip(C_grid, C_grid_ref_gas)])          
            Delta_E_grid_in.append([100*(g_in/g_ref-1) for g_in, g_ref in zip(E_grid_in, E_grid_in_ref)])


    P_grid_ref = [Case_0[0][12][str(node+1)].tolist() for node in Case_n[-1][3]]
    Qdot_boiler_ref = [Case_0[0][23][str(node+1)].tolist() for node in Case_n[-1][3]]
    P_grid = [Case_n[-1][12][str(node+1)].tolist() for node in Case_n[-1][3]]            
    
    C_grid_ref = [sum([dt*p*c for p,c in zip(profile, E_prices)]) for profile in P_grid_ref]
    C_grid_ref_gas = [sum([dt*(p*c + gas_price*q/(1000*12.7)) for p,c,q in zip(profile, E_prices, profile_gas)]) for profile,profile_gas in zip(P_grid_ref, Qdot_boiler_ref)]
    C_grid = [sum([dt*p*c for p,c in zip(profile, E_prices)])+precharge_TESS for profile in P_grid]
    
    E_grid_ref = [dt*sum(Case_0[0][12][str(node+1)].tolist()) for node in Case_n[-1][3]]
    E_grid = [dt*sum(Case_n[-1][12][str(node+1)].tolist()) for node in Case_n[-1][3]]
    
    E_grid_in_ref = [dt*sum([i for i in Case_0[0][12][str(node+1)].tolist() if i>0]) for node in Case_n[-1][3]]            
    E_grid_in = [dt*sum([i for i in Case_n[-1][12][str(node+1)].tolist() if i>0]) for node in Case_n[-1][3]]

    Delta_Cost.append([100*(c/c_ref-1) for c,c_ref in zip(C_grid, C_grid_ref)])            
    Delta_Cost_gas.append([100*(c/c_ref-1) for c,c_ref in zip(C_grid, C_grid_ref_gas)])   
    Delta_E_grid_in.append([100*(g_in/g_ref-1) for g_in, g_ref in zip(E_grid_in, E_grid_in_ref)])

            
    if Case_1:
        for penetration in range(penetrations):
            for case in range(cases):        
#                P_PV.append([Case_n[penetration*penetrations + penetration + case][7][str(node+1)].tolist() for node in Case_n[penetration*penetrations + penetration + case][3]])            
#                P_PV = [Case_n[penetration*penetrations + penetration + case][7][str(node+1)].tolist() for node in Case_n[penetration*penetrations + penetration + case][3]]
#                E_PV.append([dt*sum(Case_n[penetration*penetrations + penetration + case][7][str(node+1)].tolist()) for node in Case_n[penetration*penetrations + penetration + case][3]])        
                E_PV_ref = [dt*sum(Case_1[penetration*penetrations + penetration + case][7][str(node+1)].tolist()) for node in Case_n[penetration*penetrations + penetration + case][3]]    
                E_PV = [dt*sum(Case_n[penetration*penetrations + penetration + case][7][str(node+1)].tolist()) for node in Case_n[penetration*penetrations + penetration + case][3]]    
                Delta_E_PV_curtailed.append([100*(pv/pv_ref-1) for pv, pv_ref in zip(E_PV, E_PV_ref)])
        
                E_grid_out_ref = [dt*sum([i for i in Case_1[penetration*penetrations + penetration + case][12][str(node+1)].tolist() if i<0]) for node in Case_n[penetration*penetrations + penetration + case][3]]
                E_grid_out = [dt*sum([i for i in Case_n[penetration*penetrations + penetration + case][12][str(node+1)].tolist() if i<0]) for node in Case_n[penetration*penetrations + penetration + case][3]]
                Delta_E_grid_out.append([100*(g_out/g_ref-1) for g_out, g_ref in zip(E_grid_out, E_grid_out_ref)])


        E_PV_ref = [dt*sum(Case_1[-1][7][str(node+1)].tolist()) for node in Case_n[-1][3]]    
        E_PV = [dt*sum(Case_n[-1][7][str(node+1)].tolist()) for node in Case_n[-1][3]]        
        Delta_E_PV_curtailed.append([100*(pv/pv_ref-1) for pv, pv_ref in zip(E_PV, E_PV_ref)])
        
        E_grid_out_ref = [dt*sum([i for i in Case_1[-1][12][str(node+1)].tolist() if i<0]) for node in Case_n[-1][3]]
        E_grid_out = [dt*sum([i for i in Case_n[-1][12][str(node+1)].tolist() if i<0]) for node in Case_n[-1][3]]
        Delta_E_grid_out.append([100*(g_out/g_ref-1) for g_out, g_ref in zip(E_grid_out, E_grid_out_ref)])

    if Case_2:
        for penetration in range(penetrations):
            for case in range(cases):              
#                E_HP_ref = [dt*sum([i for i in Case_2[penetration*penetrations + penetration + case][9][str(node+1)].tolist()]) for node in Case_n[penetration*penetrations + penetration + case][3]]
#                E_HP = [dt*sum([i for i in Case_n[penetration*penetrations + penetration + case][9][str(node+1)].tolist()]) for node in Case_n[penetration*penetrations + penetration + case][3]]
#                Delta_E_HP.append([100*(hp/hp_ref-1) for hp, hp_ref in zip(E_HP, E_HP_ref)])
#
#        E_HP_ref = [dt*sum([i for i in Case_2[-1][9][str(node+1)].tolist()]) for node in Case_n[-1][3]]
#        E_HP = [dt*sum([i for i in Case_n[-1][9][str(node+1)].tolist()]) for node in Case_n[-1][3]]
#        Delta_E_HP.append([100*(hp/hp_ref-1) for hp, hp_ref in zip(E_HP, E_HP_ref)])
                E_HP_ref = [dt*sum([i+j for i,j in zip(Case_2[penetration*penetrations + penetration + case][9][str(node+1)].tolist(),Case_2[penetration*penetrations + penetration + case][10][str(node+1)].tolist())]) for node in Case_n[penetration*penetrations + penetration + case][3]]
                E_HP = [dt*sum([i+j for i,j in zip(Case_n[penetration*penetrations + penetration + case][9][str(node+1)].tolist(),Case_n[penetration*penetrations + penetration + case][10][str(node+1)].tolist())]) for node in Case_n[penetration*penetrations + penetration + case][3]]
                Delta_E_HP.append([100*(hp/hp_ref-1) for hp, hp_ref in zip(E_HP, E_HP_ref)])

        E_HP_ref = [dt*sum([i+j for i,j in zip(Case_2[-1][9][str(node+1)].tolist(),Case_2[-1][10][str(node+1)].tolist())]) for node in Case_n[-1][3]]
        E_HP = [dt*sum([i+j for i,j in zip(Case_n[-1][9][str(node+1)].tolist(),Case_n[-1][10][str(node+1)].tolist())]) for node in Case_n[-1][3]]
        Delta_E_HP.append([100*(hp/hp_ref-1) for hp, hp_ref in zip(E_HP, E_HP_ref)])                

    if Case_3:
        for penetration in range(penetrations):
            for case in range(cases):        
                E_BESS_in_ref = [dt*sum([i for i in Case_3[penetration*penetrations + penetration + case][11][str(node+1)].tolist() if i<0]) for node in Case_n[penetration*penetrations + penetration + case][3]]
                E_BESS_in = [dt*sum([i for i in Case_n[penetration*penetrations + penetration + case][11][str(node+1)].tolist() if i<0]) for node in Case_n[penetration*penetrations + penetration + case][3]]
                Delta_E_BESS_in.append([100*(BESS_in/BESS_in_ref-1) for BESS_in, BESS_in_ref in zip(E_BESS_in, E_BESS_in_ref)])

                C_BESS_ref = [Case_3[penetration*penetrations + penetration + case][16][str(node+1)].tolist()[-1] for node in Case_n[penetration*penetrations + penetration + case][3]]
                C_BESS =  [Case_n[penetration*penetrations + penetration + case][16][str(node+1)].tolist()[-1] for node in Case_n[penetration*penetrations + penetration + case][3]]
                Delta_C_BESS.append([100*(c/c_ref-1) for c, c_ref in zip(C_BESS, C_BESS_ref)])

        E_BESS_in_ref = [dt*sum([i for i in Case_3[-1][11][str(node+1)].tolist() if i<0]) for node in Case_n[-1][3]]
        E_BESS_in = [dt*sum([i for i in Case_n[-1][11][str(node+1)].tolist() if i<0]) for node in Case_n[-1][3]]
        Delta_E_BESS_in.append([100*(BESS_in/BESS_in_ref-1) for BESS_in, BESS_in_ref in zip(E_BESS_in, E_BESS_in_ref)])
    
        C_BESS_ref = [Case_3[-1][16][str(node+1)].tolist()[-1] for node in Case_n[-1][3]]
        C_BESS =  [Case_n[-1][16][str(node+1)].tolist()[-1] for node in Case_n[-1][3]]
        Delta_C_BESS.append([100*(c/c_ref-1) for c, c_ref in zip(C_BESS, C_BESS_ref)])
    
    if Case_4:
        for penetration in range(penetrations):
            for case in range(cases): 
                E_HP_TESS_ref = [dt*sum([i for i in Case_4[penetration*penetrations + penetration + case][10][str(node+1)].tolist()]) for node in Case_n[penetration*penetrations + penetration + case][3]]
                E_HP_TESS = [dt*sum([i for i in Case_n[penetration*penetrations + penetration + case][10][str(node+1)].tolist()]) for node in Case_n[penetration*penetrations + penetration + case][3]]
                Delta_E_HP_TESS.append([100*(e/e_ref-1) for e, e_ref in zip(E_HP_TESS, E_HP_TESS_ref)])

                Qdot_TESS_ref = [dt*sum([i for i in Case_4[penetration*penetrations + penetration + case][20][str(node+1)].tolist()]) for node in Case_n[penetration*penetrations + penetration + case][3]]
                Qdot_TESS = [dt*sum([i for i in Case_n[penetration*penetrations + penetration + case][20][str(node+1)].tolist()]) for node in Case_n[penetration*penetrations + penetration + case][3]]
                Delta_Qdot_TESS.append([100*(q/q_ref-1) for q, q_ref in zip(Qdot_TESS, Qdot_TESS_ref)])                


        E_HP_TESS_ref = [dt*sum([i for i in Case_4[-1][10][str(node+1)].tolist()]) for node in Case_n[-1][3]]
        E_HP_TESS = [dt*sum([i for i in Case_n[-1][10][str(node+1)].tolist()]) for node in Case_n[-1][3]]
        Delta_E_HP_TESS.append([100*(e/e_ref-1) for e, e_ref in zip(E_HP_TESS, E_HP_TESS_ref)])

        Qdot_TESS_ref = [dt*sum([i for i in Case_4[-1][20][str(node+1)].tolist()]) for node in Case_n[-1][3]]
        Qdot_TESS = [dt*sum([i for i in Case_n[-1][20][str(node+1)].tolist()]) for node in Case_n[-1][3]]
        Delta_Qdot_TESS.append([100*(q/q_ref-1) for q, q_ref in zip(Qdot_TESS, Qdot_TESS_ref)])                  
                
    
    if verbose:
        print('Delta Cost')
        print('Penetration: 20% ', min(min(Delta_Cost[0:5])), max(max(Delta_Cost[0:5])), len(Delta_Cost[0:5]))
        print('Penetration: 40% ', min(min(Delta_Cost[5:10])), max(max(Delta_Cost[5:10])), len(Delta_Cost[5:10]))
        print('Penetration: 60% ', min(min(Delta_Cost[10:15])), max(max(Delta_Cost[10:15])), len(Delta_Cost[10:15]))
        print('Penetration: 80% ', min(min(Delta_Cost[15:20])), max(max(Delta_Cost[15:20])), len(Delta_Cost[15:20]))

        print('Delta Cost with Gas')
        print('Penetration: 20% ', min(min(Delta_Cost_gas[0:5])), max(max(Delta_Cost_gas[0:5])), len(Delta_Cost_gas[0:5]))
        print('Penetration: 40% ', min(min(Delta_Cost_gas[5:10])), max(max(Delta_Cost_gas[5:10])), len(Delta_Cost_gas[5:10]))
        print('Penetration: 60% ', min(min(Delta_Cost_gas[10:15])), max(max(Delta_Cost_gas[10:15])), len(Delta_Cost_gas[10:15]))
        print('Penetration: 80% ', min(min(Delta_Cost_gas[15:20])), max(max(Delta_Cost_gas[15:20])), len(Delta_Cost_gas[15:20]))        
  
        print('Delta E grid in')
        print('Penetration: 20% ', min(min(Delta_E_grid_in[0:5])), max(max(Delta_E_grid_in[0:5])), len(Delta_E_grid_in[0:5]))
        print('Penetration: 40% ', min(min(Delta_E_grid_in[5:10])), max(max(Delta_E_grid_in[5:10])), len(Delta_E_grid_in[5:10]))
        print('Penetration: 60% ', min(min(Delta_E_grid_in[10:15])), max(max(Delta_E_grid_in[10:15])), len(Delta_E_grid_in[10:15]))
        print('Penetration: 80% ', min(min(Delta_E_grid_in[15:20])), max(max(Delta_E_grid_in[15:20])), len(Delta_E_grid_in[15:20]))
 
    if verbose and Case_1:      
        print('Delta E grid out')
        print('Penetration: 20% ', min(min(Delta_E_grid_out[0:5])), max(max(Delta_E_grid_out[0:5])), len(Delta_E_grid_out[0:5]))
        print('Penetration: 40% ', min(min(Delta_E_grid_out[5:10])), max(max(Delta_E_grid_out[5:10])), len(Delta_E_grid_out[5:10]))
        print('Penetration: 60% ', min(min(Delta_E_grid_out[10:15])), max(max(Delta_E_grid_out[10:15])), len(Delta_E_grid_out[10:15]))
        print('Penetration: 80% ', min(min(Delta_E_grid_out[15:20])), max(max(Delta_E_grid_out[15:20])), len(Delta_E_grid_out[15:20]))

        print('E PV')
        print('Penetration: 20% ', min(min(Delta_E_PV_curtailed[0:5])), max(max(Delta_E_PV_curtailed[0:5])), len(Delta_E_PV_curtailed[0:5]))
        print('Penetration: 40% ', min(min(Delta_E_PV_curtailed[5:10])), max(max(Delta_E_PV_curtailed[5:10])), len(Delta_E_PV_curtailed[5:10]))
        print('Penetration: 60% ', min(min(Delta_E_PV_curtailed[10:15])), max(max(Delta_E_PV_curtailed[10:15])), len(Delta_E_PV_curtailed[10:15]))
        print('Penetration: 80% ', min(min(Delta_E_PV_curtailed[15:20])), max(max(Delta_E_PV_curtailed[15:20])), len(Delta_E_PV_curtailed[15:20]))

    if verbose and Case_2:      
        print('Delta E HP')
        print('Penetration: 20% ', min(min(Delta_E_HP[0:5])), max(max(Delta_E_HP[0:5])), len(Delta_E_HP[0:5]))
        print('Penetration: 40% ', min(min(Delta_E_HP[5:10])), max(max(Delta_E_HP[5:10])), len(Delta_E_HP[5:10]))
        print('Penetration: 60% ', min(min(Delta_E_HP[10:15])), max(max(Delta_E_HP[10:15])), len(Delta_E_HP[10:15]))
        print('Penetration: 80% ', min(min(Delta_E_HP[15:20])), max(max(Delta_E_HP[15:20])), len(Delta_E_HP[15:20]))

    if verbose and Case_3:      
        print('Delta E_BESS_in')
        print('Penetration: 20% ', min(min(Delta_E_BESS_in[0:5])), max(max(Delta_E_BESS_in[0:5])), len(Delta_E_BESS_in[0:5]))
        print('Penetration: 40% ', min(min(Delta_E_BESS_in[5:10])), max(max(Delta_E_BESS_in[5:10])), len(Delta_E_BESS_in[5:10]))
        print('Penetration: 60% ', min(min(Delta_E_BESS_in[10:15])), max(max(Delta_E_BESS_in[10:15])), len(Delta_E_BESS_in[10:15]))
        print('Penetration: 80% ', min(min(Delta_E_BESS_in[15:20])), max(max(Delta_E_BESS_in[15:20])), len(Delta_E_BESS_in[15:20]))

        print('Delta C_BESS')
        print('Penetration: 20% ', min(min(Delta_C_BESS[0:5])), max(max(Delta_C_BESS[0:5])), len(Delta_C_BESS[0:5]))
        print('Penetration: 40% ', min(min(Delta_C_BESS[5:10])), max(max(Delta_C_BESS[5:10])), len(Delta_C_BESS[5:10]))
        print('Penetration: 60% ', min(min(Delta_C_BESS[10:15])), max(max(Delta_C_BESS[10:15])), len(Delta_C_BESS[10:15]))
        print('Penetration: 80% ', min(min(Delta_C_BESS[15:20])), max(max(Delta_C_BESS[15:20])), len(Delta_C_BESS[15:20]))

    if verbose and Case_4:
        print('Delta E_HP_TESS_ref')
        print('Penetration: 20% ', min(min(Delta_E_HP_TESS[0:5])), max(max(Delta_E_HP_TESS[0:5])), len(Delta_E_HP_TESS[0:5]))
        print('Penetration: 40% ', min(min(Delta_E_HP_TESS[5:10])), max(max(Delta_E_HP_TESS[5:10])), len(Delta_E_HP_TESS[5:10]))
        print('Penetration: 60% ', min(min(Delta_E_HP_TESS[10:15])), max(max(Delta_E_HP_TESS[10:15])), len(Delta_E_HP_TESS[10:15]))
        print('Penetration: 80% ', min(min(Delta_E_HP_TESS[15:20])), max(max(Delta_E_HP_TESS[15:20])), len(Delta_E_HP_TESS[15:20]))

        print('Delta Qdot_TESS')
        print('Penetration: 20% ', min(min(Delta_Qdot_TESS[0:5])), max(max(Delta_Qdot_TESS[0:5])), len(Delta_Qdot_TESS[0:5]))
        print('Penetration: 40% ', min(min(Delta_Qdot_TESS[5:10])), max(max(Delta_Qdot_TESS[5:10])), len(Delta_Qdot_TESS[5:10]))
        print('Penetration: 60% ', min(min(Delta_Qdot_TESS[10:15])), max(max(Delta_Qdot_TESS[10:15])), len(Delta_Qdot_TESS[10:15]))
        print('Penetration: 80% ', min(min(Delta_Qdot_TESS[15:20])), max(max(Delta_Qdot_TESS[15:20])), len(Delta_Qdot_TESS[15:20]))
        
        

    if plots:
        from numpy import arange
        import matplotlib.pylab as plt

        plt.rcParams.update({
        #    "text.usetex": True,
            "font.family": "Times New Roman",
            'font.size': 16
        })
    
        E_load_ref = [i/1000 for i in E_load_ref]
        
        E_grid_out = [-dt*sum([i for i in Case_n[-1][12][str(node+1)].tolist() if i<0]) for node in Case_n[-1][3]]

#        plt.figure(constrained_layout=True)
#        plt.scatter(E_load_ref, [i*7/365 for i in E_load_ref])
#        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
#        plt.xlim((0, max(E_load_ref)+1))
#        plt.title('Weekly consumption vs. Yearly consumption')
#        plt.xlabel('Node yearly consumption [MWh]')
#        plt.ylabel('Node weekly consumption [MWh]')
#        plt.show() 
#    
#        plt.figure(constrained_layout=True)
#        plt.scatter(E_load_ref, C_grid)
#        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
#        plt.xlim((0, max(E_load_ref)+1))
#        plt.title('Cost vs. Consumption')
#        plt.xlabel('Node yearly consumption [MWh]')
#        plt.ylabel('Cost')
#        plt.show()              

#        plt.figure(constrained_layout=True)
#        plt.grid()                
#        plt.scatter(E_load_ref, Delta_Cost_gas[-1], color = 'r', label = 'Electricity + gas')
#        plt.scatter(E_load_ref, Delta_Cost[-1], label = 'Electricity')
#        plt.legend(loc='upper right')
#        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
#        plt.xlim((0, max(E_load_ref)+1))
#        plt.title('Change in cost vs. Consumption')
#        plt.xlabel('Node yearly consumption [MWh]')
#        plt.ylabel('Change in cost [%]')
#        plt.show()    
        
        fig, ax1 = plt.subplots()
        ax1 = plt.gca()
        ax1.scatter(E_load_ref, Delta_Cost[-1], label = 'Electricity')
        ax1.set_xlabel('Node yearly consumption [MWh]')       
        ax1.set_xlim((0, max(E_load_ref)+1))
        ax1.set_xticks([i for i in range(16)])
        ax1.tick_params(axis='y', labelcolor='tab:blue')
        ax1.set_ylabel('Change in electricity cost [%]')        
        ax1.grid()
        
        ax2 = plt.gca()
        ax2 = ax1.twinx()
        ax2.scatter(E_load_ref, Delta_Cost_gas[-1], color = 'r', label = 'Electricity + gas')
        ax2.set_ylabel('Change in energy cost [%]')        
        ax2.tick_params(axis='y', labelcolor='r')
#        ax2.grid(None)
        
#        plt.xticks(E_load_ref, E_load_ref)
#        plt.grid(True)
        fig.tight_layout()  
        plt.show()

        plt.figure(constrained_layout=True)
        plt.grid()        
        plt.scatter(E_load_ref, Delta_E_grid_in[-1])
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('Change in energy consumption vs. consumption')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('Change in energy consumption [%]')
        plt.show() 

#        plt.figure(constrained_layout=True)
#        plt.grid()        
#        plt.scatter(E_load_ref, [(e_out/1000)/(l*7/365) for e_out,l in zip(E_grid_out, E_load_ref)])
#        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
#        plt.xlim((0, max(E_load_ref)+1))
#        plt.title('Ratio of energy injection vs. consumption')
#        plt.xlabel('Node yearly consumption [MWh]')
#        plt.ylabel('Ratio of energy injection')
#        plt.show()        

#        plt.figure(constrained_layout=True)
#        plt.grid()        
#        plt.scatter(E_load_ref, E_PV)
#        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
#        plt.xlim((0, max(E_load_ref)+1))
#        plt.title('')
#        plt.xlabel('Node yearly consumption [MWh]')
#        plt.ylabel('')
#        plt.show()    
        
        if Case_1:
            E_PV = [i/1000 for i in E_PV] 
#            E_HP = [dt*sum([i for i in Case_2[-1][9][str(node+1)].tolist()]) for node in Case_n[-1][3]]
            E_HP = [dt*sum([(i+j)/1000 for i,j in zip(Case_n[-1][9][str(node+1)].tolist(),Case_n[-1][10][str(node+1)].tolist())]) for node in Case_n[-1][3]]
            
#            plt.figure(constrained_layout=True)
#            plt.grid()            
#            plt.scatter(E_load_ref, E_PV)
#            plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
#            plt.xlim((0, max(E_load_ref)+1))
#            plt.title('Generation vs. Consumption')
#            plt.xlabel('Node yearly consumption [MWh]')
#            plt.ylabel('PV generation [MWh]')

#            plt.show()   

#            plt.figure(constrained_layout=True)
#            plt.grid()            
#            plt.scatter(E_load_ref, [pv/(l*7/365) for pv,l in zip(E_PV, E_load_ref)])
#            plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
#            plt.xlim((0, max(E_load_ref)+1))
#            plt.title('Generation-consumption ratio vs. consumption')
#            plt.xlabel('Node yearly consumption [MWh]')
#            plt.ylabel('Generation-consumption ratio')

#            plt.show()       

            plt.figure(constrained_layout=True)
            plt.grid()            
            plt.scatter(E_load_ref, Delta_E_grid_out[-1])
            plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
            plt.xlim((0, max(E_load_ref)+1))
            plt.title('Change in energy injection vs. consumption')
            plt.xlabel('Node yearly consumption [MWh]')
            plt.ylabel('Change in energy injection [%]')

            plt.show()

            plt.figure(constrained_layout=True)
            plt.grid()            
#            plt.scatter(E_load_ref, [Delta_E_HP])
            plt.scatter(E_load_ref, [hp/(l*7/365) for hp,l in zip(E_HP, E_load_ref)])
            plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
            plt.xlim((0, max(E_load_ref)+1))
            plt.title('Ratio of HP consumption vs. Consumption')
            plt.xlabel('Node yearly consumption [MWh]')
            plt.ylabel('Ratio of HP consumption')

#            plt.show()   
            

        if Case_2:            
#            Delta_E_HP = [i/1000 for i in Delta_E_HP[-1]] 
            E_BESS_in = [-dt*sum([i for i in Case_n[-1][11][str(node+1)].tolist() if i<0]) for node in Case_n[-1][3]]
            E_BESS_out = [dt*sum([i for i in Case_n[-1][11][str(node+1)].tolist() if i>0]) for node in Case_n[-1][3]]
            
            plt.figure(constrained_layout=True)
            plt.grid()            
            plt.scatter(E_load_ref, Delta_E_HP[-1])
            plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
            plt.xlim((0, max(E_load_ref)+1))
            plt.title('Change in HP consumption vs. Consumption')
            plt.xlabel('Node yearly consumption [MWh]')
            plt.ylabel('Change in HP consumption [%]')
            plt.show()  

            plt.figure(constrained_layout=True)
            plt.grid()            
#            plt.scatter(E_load_ref, [pv/(l*7/365) for pv,l in zip(E_PV, E_load_ref)])
#            plt.scatter(E_load_ref, E_BESS_in, c='red', label = 'Charge')
            plt.scatter(E_load_ref, E_BESS_out) # , c='blue', label = 'Discharge'
            plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
            plt.xlim((0, max(E_load_ref)+1))
            plt.title('Energy supplied by the BESS vs. Consumption')
            plt.xlabel('Node yearly consumption [MWh]')
            plt.ylabel('Energy supply by the BESS [kWh]')
#            plt.legend(loc='lower center', ncol=1)
            plt.show() 

#        if Case_3:

        if Case_4:
            plt.figure(constrained_layout=True)
            plt.grid()            
            plt.scatter(E_load_ref, Delta_E_HP_TESS[-1])
            plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
            plt.xlim((0, max(E_load_ref)+1))
            plt.title('Change in energy used to charge the TESS vs. consumption')
            plt.xlabel('Node yearly consumption [MWh]')
            plt.ylabel('Change in energy used to charge the TESS [%]')
            plt.show()

            plt.figure(constrained_layout=True)
            plt.grid()            
            plt.scatter(E_load_ref, Delta_Qdot_TESS[-1])
            plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
            plt.xlim((0, max(E_load_ref)+1))
            plt.title('Change in heat delivered by the TESS vs. consumption')
            plt.xlabel('Node yearly consumption [MWh]')
            plt.ylabel('Change in heat delivered by the TESS [%]')
            plt.show() 

    if Case_6:
        T_0 = []
        T_1 = []
        T_2 = []
        T_3 = []
        T_4 = []
        T_5 = []
        T_6 = []
        T_7 = []

        for node in Case_n[-1][3]:
            T_0 += Case_0[-1][25][str(node+1)].tolist()
            T_1 += Case_1[-1][25][str(node+1)].tolist()
            T_2 += Case_2[-1][25][str(node+1)].tolist()
            T_3 += Case_3[-1][25][str(node+1)].tolist()
            T_5 += Case_4[-1][25][str(node+1)].tolist()
            T_6 += Case_5[-1][25][str(node+1)].tolist()
            T_7 += Case_6[-1][25][str(node+1)].tolist()
            T_4 += Case_n[-1][25][str(node+1)].tolist()

        T_0 = [i-273 for i in T_0]
        T_1 = [i-273 for i in T_1]
        T_2 = [i-273 for i in T_2]
        T_3 = [i-273 for i in T_3]
        T_4 = [i-273 for i in T_4]
        T_5 = [i-273 for i in T_5]
        T_6 = [i-273 for i in T_6]
        T_7 = [i-273 for i in T_7]
            
        plt.figure(constrained_layout=True)
#        plt.grid()          
        plt.boxplot([T_0, T_1, T_2, T_3, T_4, T_5, T_6, T_7], labels=['0', '1', '2', '3', '4', '5', '6', '7'], showfliers=False)
        plt.title('Temperature distribution')
        plt.xlabel('Cases')
        plt.ylabel('Indoor temperature, T$_{in}$, [°C]')
        plt.show() 
        
#        from numpy import linspace
#        
#        bins = linspace(16, 20, 40)
#
#        plt.figure(constrained_layout=True)
#        plt.hist([T_0, T_1, T_2, T_3, T_4, T_5, T_6], bins, label=['0', '1', '2', '3', '4', '5', '6'])           
##        plt.hist(T_0, bins, alpha=0.5, label='Case 0')
##        plt.hist(T_1, bins, alpha=0.5, label='Case 1')
##        plt.hist(T_2, bins, alpha=0.5, label='Case 2')
##        plt.hist(T_3, bins, alpha=0.5, label='Case 3')
##        plt.hist(T_4, bins, alpha=0.5, label='Case 4')
##        plt.hist(T_5, bins, alpha=0.5, label='Case 5')
##        plt.hist(T_6, bins, alpha=0.5, label='Case 6')
#        plt.xlabel('Indoor temperature, T$_{in}$, [°C]')
#        plt.ylabel('Ocurrences')
#        plt.legend(loc='upper left')
#        plt.show()        
                
            
#    return [E_load_ref, E_PV]# [Delta_Cost, Delta_E_grid_in, Delta_E_grid_out, Delta_P_PV_curtailed]

###############################################################################

def TESS_SoC(T_TESS, T_min = 50+273, T_max = 90+273):
    
    return (100/(T_max - T_min))*T_TESS + (-(100/(T_max - T_min))*T_min)

###############################################################################

def Plot_objectives(case_n, selected_nodes_list, penetration = 20, case = 0, start = 0, end = 24*4*7, i = -1, penetrations = [10*i for i in arange(1,10,1)], Network_headers = ['Van.Nummer', 'Naar.Nummer', 'Lengte', 'Kabeltype', 'GM', 'House_Type'], dt = 0.25):
    from numpy import array, linspace
    import matplotlib.pylab as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })      

    if i != -1:
        if penetration != 100:
            i = 5*(penetrations.index(penetration)) + case
        else:
            i = -1
        
    penetration = penetrations.index(penetration)
    
    # Temperature variation from setpoint
    
    T_set_day = [273 + 20 - 3]*int((6-0)*4) + [273 + 20]*int((22-6)*4) + [273 + 20 - 3]*int((24-22)*4)
    T_set = array(T_set_day*int(len(range(end))/24/4))
    
    T_incompliance_studio = []
    T_incompliance_apartment = []
    T_incompliance_stand_alone = []

    # Amount of energy consumed from the grid and its cost
    E_grid_studio = []
    E_grid_apartment = []
    E_grid_stand_alone = []
    c_grid_studio = []
    c_grid_apartment = []
    c_grid_stand_alone = []    
    
    
    for node in selected_nodes_list[penetration][case]:
        if case_n[i][30][Network_headers[5]][node+1] == 'studio':
            # print('Node: ', node, ', studio, ', case_n[i][8][str(node+1)].sum()*dt*365/7)
            T_incompliance_studio += [T-s for T,s in zip(case_n[i][25][str(node+1)].tolist(), T_set)]
            E_grid_studio.append(case_n[i][12][str(node+1)].sum()*dt)
            c_grid_studio.append(dt*sum([P*c for P,c in zip(case_n[i][12][str(node+1)], case_n[i][29]['Current'])]))
            
        elif case_n[i][30][Network_headers[5]][node+1] == 'apartment':
            # print('Node: ', node, ', apartment', case_n[i][8][str(node+1)].sum()*dt*365/7)
            T_incompliance_apartment += [T-s for T,s in zip(case_n[i][25][str(node+1)].tolist(), T_set)]
            E_grid_apartment.append(case_n[i][12][str(node+1)].sum()*dt)
            c_grid_apartment.append(dt*sum([P*c for P,c in zip(case_n[i][12][str(node+1)], case_n[i][29]['Current'])]))
            
        elif case_n[i][30][Network_headers[5]][node+1] == 'stand_alone':
            # print('Node: ', node, ', stand_alone', case_n[i][8][str(node+1)].sum()*dt*365/7)
            T_incompliance_stand_alone += [T-s for T,s in zip(case_n[i][25][str(node+1)].tolist(), T_set)]
            E_grid_stand_alone.append(case_n[i][12][str(node+1)].sum()*dt)
            c_grid_stand_alone.append(dt*sum([P*c for P,c in zip(case_n[i][12][str(node+1)], case_n[i][29]['Current'])]))
    
    
    bins = linspace(-5, 5, 10)
    plt.figure(constrained_layout = True)
    # plt.hist(T_incompliance_studio, bins, alpha=0.5, label='Studio', color = 'b', density = True)
    # plt.hist(T_incompliance_apartment, bins, alpha=0.5, label='Apartment', color = 'r', density = True)
    # plt.hist(T_incompliance_stand_alone, bins, alpha=0.5, label='Stand alone', color = 'k', density = True)
    plt.hist([T_incompliance_studio, T_incompliance_apartment, T_incompliance_stand_alone], bins, alpha=0.5, label=['Studio', 'Apartment', 'Stand alone'], color = ['b', 'r', 'k'], density = True)
    plt.legend(loc='upper left')
    plt.show()

    
    # Amount of energy consumed from the grid and its cost
    # bins = linspace(0, 1000, 10)
    plt.figure(constrained_layout = True)
    # plt.hist(E_grid_studio, alpha=0.5, label='Studio', color = 'b', density = True)
    # plt.hist(E_grid_apartment, alpha=0.5, label='Apartment', color = 'r', density = True)
    # plt.hist(E_grid_stand_alone, alpha=0.5, label='Stand alone', color = 'k', density = True)
    plt.hist([E_grid_studio, E_grid_apartment, E_grid_stand_alone], alpha=0.5, label=['Studio', 'Apartment', 'Stand alone'], color = ['b', 'r', 'k'])
    plt.legend(loc='upper left')
    plt.show()

    # bins = linspace(0, 1000, 10)
    plt.figure(constrained_layout = True)
    # plt.hist(E_grid_studio, alpha=0.5, label='Studio', color = 'b', density = True)
    # plt.hist(E_grid_apartment, alpha=0.5, label='Apartment', color = 'r', density = True)
    # plt.hist(E_grid_stand_alone, alpha=0.5, label='Stand alone', color = 'k', density = True)
    plt.hist([c_grid_studio, c_grid_apartment, c_grid_stand_alone], alpha=0.5, label=['Studio', 'Apartment', 'Stand alone'], color = ['b', 'r', 'k'])
    plt.legend(loc='upper left')
    plt.show()    
    
    
    # Variation from the DSO setpoint and overall voltage distribution

###############################################################################
            
def Plot_house_behaviour(case_n, node, penetration = 20, case = 0, i = -1, start = 0, end = 24*4*7, penetrations = [10*i for i in arange(1,10,1)]):
    from numpy import array#linspace, sin, pi
    import matplotlib.pylab as plt

    plt.rcParams.update({
    #    "text.usetex": True,
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
#    T_amb =
    SoC_TESS = [TESS_SoC(i) for i in case_n[i][21][str(node)].tolist()]
    
    
    fig, axs = plt.subplots(7, 2, sharex=True) # , constrained_layout=True

    # Electric
    axs[0, 0].grid()
    axs[0, 0].plot(P_L) # Load
    axs[0, 0].set_ylabel('$P_{L}$, [kW]')
    axs[1, 0].grid()
    axs[1, 0].plot(P_PV_av, color='r', label='av') # PV av
    axs[1, 0].plot(P_PV, color='b', label='PV') # PV av, PV
    axs[1, 0].set_ylabel('$P_{PV}$, [kW]')
    axs[1, 0].legend(loc='upper left', ncol = 2, fontsize=10)
    axs[2, 0].grid()
    axs[2, 0].plot(P_HP_TESS, color='b', label='TESS') # TESS    
    axs[2, 0].plot(P_HP, color='r', label='L') # HP    
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
#    axs[6, 0].set_xticks([i*96 for i in range(8)])
#    axs[6, 0].set_xticklabels(['01/02', '02/02', '03/02', '04/02', '05/02', '06/02', '07/02', '08/02'])
#    axs[6, 0].set_xticklabels(['26/06', '27/06', '28/06', '29/06', '30/06', '01/07', '02/07', '03/07']) 
#    axs[6, 0].set_xticklabels(['01/02', '02/02', '03/02', '04/02', '05/02', '06/02', '07/02', '08/02'])
#    axs[6, 0].set_xticklabels(['19/07', '20/07', '21/07', '22/07', '23/07', '24/07', '25/07', '26/07'])
    
    # Thermal
    axs[0, 1].grid()
    axs[0, 1].plot(Qdot_L)
    axs[0, 1].set_ylabel('$\dot{Q}_{L}$, [kW]')    
    axs[1, 1].grid()
    axs[1, 1].plot(Qdot_Boiler)
    axs[1, 1].set_ylabel('$\dot{Q}_{B}$, [kW]')    
    axs[2, 1].grid()    
    axs[2, 1].plot(Qdot_HP_TESS, color='b', label='TESS') # TESS    
    axs[2, 1].plot(Qdot_HP, color='r', label='L') # HP
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
#    axs[5, 1].plot(y, color='y', label='amb')
    axs[5, 1].set_ylabel('T, [°C]')
    axs[5, 1].legend(loc='lower right', ncol = 1, fontsize=10)      
    axs[6, 1].grid()
    axs[6, 1].plot(Prices)    
    axs[6, 1].set_ylabel('$\lambda$, [€/kWh]')
    axs[6, 1].set_xlim([start,end])
#    axs[6, 1].set_xticks([i*96 for i in range(8)])
#    axs[6, 1].set_xticklabels(['01/02', '02/02', '03/02', '04/02', '05/02', '06/02', '07/02', '08/02'])  
#    axs[6, 1].set_xticklabels(['01/02', '02/02', '03/02', '04/02', '05/02', '06/02', '07/02', '08/02'])
    
###############################################################################

###############################################################################
            
def Plot_grid_behaviour(case_n, case_centralized, node = 269, i = -1, start = 0, end = 24*4*7, penetrations = [10*i for i in arange(1,10,1)]):
    from numpy import array#linspace, sin, pi
    import matplotlib.pylab as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 13
    })      

#    if i != -1:
#        if penetration != 100:
#            i = 5*(penetrations.index(penetration)) + case
#        else:
#            i = -1
#    
#    P_L = case_n[i][12][str(node)].tolist()
#
#    P_BESS = case_n[i][11][str(node)].tolist()
#    SoC_BESS = [100*i for i in case_n[i][15][str(node)].tolist()]
#    
#    P_Grid = case_n[i][12][str(node)].tolist()
#    P_Grid_DSO = case_n[i][13][str(node)].tolist()

    
#    l_min_centralized = [min(df_V_n.loc[i])/400 for i in range(end)]
#    l_max_centralized = [max(df_V_n.loc[i])/400 for i in range(end)]
#
#
#    l_min = [min(Case_2_PV_HP_W[-1][27].loc[i]) for i in range(end)]
#    l_max = [max(Case_2_PV_HP_W[-1][27].loc[i]) for i in range(end)]    
        
    fig, axs = plt.subplots(3, sharex=True) # , constrained_layout=True

    axs[0].grid()
    axs[0].set_ylabel('$V_{n}$, [p.u.]')
    axs[0].plot([min(case_n[i][27].loc[ts])/400 for ts in range(end)], color='b', label='Case 2')
    axs[0].plot([max(case_n[i][27].loc[ts])/400 for ts in range(end)], color='b')
    axs[0].plot([min(case_centralized[i][15].loc[ts])/400 for ts in range(end)], color='r', label='Case 8')
    axs[0].plot([max(case_centralized[i][15].loc[ts])/400 for ts in range(end)], color='r')
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
    axs[2].set_xlim([start,end])
#    axs[2].set_xticks([i*96 for i in range(8)])
#    axs[2].set_xticklabels(['01/02', '02/02', '03/02', '04/02', '05/02', '06/02', '07/02', '08/02'])
#    axs[2].set_xticklabels(['26/06', '27/06', '28/06', '29/06', '30/06', '01/07', '02/07', '03/07']) 
#    axs[2].set_xticklabels(['01/02', '02/02', '03/02', '04/02', '05/02', '06/02', '07/02', '08/02'])
#    axs[2].set_xticklabels(['19/07', '20/07', '21/07', '22/07', '23/07', '24/07', '25/07', '26/07'])

###############################################################################


#
#def Cost_function(d_P, df_P_Grid = False, ts = False):
##    P_Grid_ts = df_P_Grid.loc[ts].tolist()
#    
#    return sum([(P)**2 for P in d_P])
##    return sum([(P/Sn)**2 for Sn,P in zip(P_Grid_ts, d_P)])
    
    

###############################################################################

def Voltage_Constraints(A, Z, df_P_Grid, d_P, ts, I_0, V_0):
    # [V, I] = Voltage_Constraints(A, Z, df_S_n_0, d_P, ts, I0.loc[ts].tolist(), V0.loc[ts].tolist())
    P_Grid_ts = df_P_Grid.loc[ts].tolist()
#    P_Grid_new = [(Sn + P)*1000 for Sn,P in zip(P_Grid_ts, d_P)]
    
    P_Grid_new = [Sn*(1 + P)*1000 for Sn,P in zip(P_Grid_ts, d_P)]
    
    return Estimate_Node_voltage_2(A, Z, P_Grid_new, I_0, V_0)

###############################################################################

def Estimate_Currents_2(A, S_n, I_0, V_0, dI = 0.01, show_plot = False):
    from numpy import array, divide, matmul, transpose, where, isnan
    
#    from pandas import isna

#    print('----------------------------------')
#    print(type(A))
#    print(type(S_n))    
    
#    if [] != where(isnan(divide(S_n, V_0))):
##    if [] != where(isna(divide(S_n, V_0))):        
#        I_n = divide(S_n, V_0)
#    else:
#        print('zeros')
#        I_n = array([S/V if V!=0 else 0 for S,V in zip(S_n, V_0)])

    I_n = divide(S_n, V_0)

    I = matmul(transpose(A),transpose(I_0)) + I_n
#    I = [a+b for a,b in zip(matmul(transpose(A),transpose(I_0)), I_n)]
#    I = I.tolist()[1:]
    
#    print('----------------------------------')
#    print(I)
#    print(type(I_0))
#    print(I_0)
    
#    i = 0
##    while max(abs(I - I_0)) > dI:
#    for i in range(10):
#        
#        I_0 = I
#        I = matmul(transpose(A),I_0) + I_n
#        
##        i+=1
##        if i>100:
##            print('Error in current estimation')
##            break
#
#    if show_plot:
#        plot_current(I)
        
    return I

###############################################################################

def Estimate_Node_voltage_2(A, Z, S_n, I_0, V_0, V_feeder = 400, dV = 0.001, dI = 0.01, show_plot = False): # V_feeder = 400
    from numpy import matmul, zeros
    
    I_0 = Estimate_Currents_2(A, S_n, I_0, V_0, dI)
    
    B = zeros(len(I_0))
    B[0] = 1

    V_n = matmul(A,V_0) - matmul(Z,I_0) + B*V_feeder
    
#    V_n = zeros(len(I_0)) + V_feeder
    
#    print(I_0)
#    print(V_n)
    
#    i = 0
#    while max(abs(V_n - V_0)) > dV:
#        
#        V_0 = V_n
#        I_0 = Estimate_Currents_2(A, S_n, I_0, V_0, dI)
#        V_n = matmul(A,V_0) - matmul(Z,I_0) + B*V_feeder
#        
#        i+=1
#        if i>1000:
#            print('Error in voltage estimation: ', max(abs(V_n - V_0)))
#            break    
#    
#    
#    if show_plot:
#      plot_voltage(V_n, V_feeder)
     
#    print(I_0)
        
    return [V_n, I_0]



###############################################################################
##############################   Example Gekko   ##############################

#def Create_setpoints():
#    from gekko import GEKKO
#    import numpy as np
#    import matplotlib.pyplot as plt 
#    
#    
#    DF_Network = pd.read_excel('Stedin_network_3.xlsx') # ['Gaia_network_2.xlsx', 'Stedin_network_2.xlsx', 'Stedin_network_3']
#    Network_headers = ['Van.Nummer', 'Naar.Nummer', 'Lengte', 'Kabeltype', 'GM', 'House_Type']
#    Wires_headers = ['Wire', 'R', 'X'] 
#    wire_file = 'Gaia_cables.xlsx'
#    csv_name = 'S_n_Stedin_301.csv'
#    
#    [A, Z] = Create_network(DF_Network, Network_headers, Wires_headers, wire_file = 'Gaia_cables.xlsx')
#    
#    
#    ts = 0
#    
#    
#    df_S_n_0 = Case_0_base_W[0][12]
#    V0 = Case_0_base_W[0][27].loc[ts].tolist()
#    I0 = Case_0_base_W[0][28].loc[ts].tolist()
#    
#    df_S_n_4 = Case_4_MCES_NA_S[-1][12]
#    V4 = Case_4_MCES_NA_S[-1][27].loc[ts].tolist()
#    I4 = Case_4_MCES_NA_S[-1][28].loc[ts].tolist()
#
#    
#    
#    # Create the Gekko model
#    m = GEKKO()    
#    
#    # Create a continuous variable    
#    delta_P_grid = m.Array(m.Var, len(df_S_n_0.loc[ts].tolist()), lb=-2, ub=2)
#    
#    
#    #[V, I] = Voltage_Constraints(A, Z, df_S_n_0, d_P, ts, I0.loc[ts].tolist(), V0.loc[ts].tolist())
#    
#    # create the problem
#    m.Minimize(Cost_function(delta_P_grid))
#     
#    
#    # Create constrain
#    m.Equations([V_n > 0.95 for V_n in Voltage_Constraints(A, Z, df_S_n_0, delta_P_grid, ts, I0, V0)[0]])
#    
#    
#    ############################   Optimize  ##################################
#    
#    m.options.SOLVER = 1 # 1=APOPT, 2=BPOPT, 3=IPOPT
#    m.solve(disp = True)


###############################################################################
##############################   Example SciPy   ##############################


# Define the objective function to minimize (profit function)
def Cost_function_1(x, *args): #, *args
    from numpy import divide, matmul, transpose
    Sn0 = args[5]    
#    return sum([(P)**2 for P in x])
    return sum([(P - Sn)**2 for P, Sn in zip(x, Sn0)])  
    
def Cost_function_2(x, *args): #, *args
    from numpy import divide, matmul, transpose
#    df_S_n_0 = args
    A = args[0]
    B = args[1]
    Z = args[2]
    V0 = args[3]
    I0 = args[4]
#    return sum([(P)**2 for P in x])
#    return sum([abs(matmul(A[i,:],V0) - matmul(Z[i,:],x) + B[i]*400 - 400) for i in range(301)])
#    return sum([(P - I)**2 for P, I in zip(x, I0)])
#    return sum([(P - Sn)**2 for P, Sn in zip(x, Sn0)])    
    return sum([abs(abs(matmul(A[i,:],V0) - matmul(Z[i,:],(matmul(transpose(A),transpose(I0)) + divide(x, [i/1000 for i in V0]))) + B[i]*400) - 400)**2 for i in range(301)])

def Cost_function(x, *args):
    
    weights = args[6]
    
    return weights[0]*Cost_function_1(x, *args) + weights[1]*Cost_function_2(x, *args)

def Find_power_setpoints(A, Z, B, V0, I0, Sn0, all_nodes, selected_nodes, pu_lim = 0.05, v_feeder = 400, weights = [1162690, 13], plotting = False):
    #import pandas as pd
    from scipy.optimize import Bounds, LinearConstraint, minimize
    #import numpy as np
    from numpy import divide, matmul, transpose, zeros
    # https://www.youtube.com/watch?v=X0LvnxSqfNk&ab_channel=KodyPowell
    
    # Define the bounds for each decision variable
#    bounds = [(-10, 10) for i in I0]
    
    # Initial guess for the decision variables
    initial_guess = [i for i in Sn0]
    
    # Define the inequality constraint functions    

    cons = []#[{'type': 'ineq', 'fun': lambda x: abs(Estimate_Node_voltage(A, Z, x, I0, V0)[0][270]) - pu_lim*400}]
#    c0_low = {'type': 'ineq', 'fun': lambda x : (-0.00083106*sum(x)**2)/(114**2) - 0.03378*sum(x)/114 + 0.9997 - (1-pu_lim)}
#    cons.append(c0_low)
#    c0_high = {'type': 'ineq', 'fun': lambda x : (0.00083106*sum(x)**2)/(114**2) + 0.03378*sum(x)/114 - 0.9997 + (1+pu_lim)}
#    cons.append(c0_high)
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
    
#    start_time = time.time()    # The timer is initializad.
    opt = minimize(Cost_function_1, initial_guess, method = 'SLSQP', args = (A, B, Z, V0, I0, Sn0, weights), constraints = cons, options={'disp': True}) #, method = 'SLSQP', bounds = bounds, constraints = cons, args = (A, Z, df_S_n_0, ts, I0, V0), options={'maxiter': 1000}
#    end_time = time.time()    # The timer is finished
    
#    print('Computation time: ', (end_time - start_time)/60)
    
    # Print the result
#    print(opt.message)
#    print(opt.success)
#    print(cons)

    if plotting:
        import matplotlib.pyplot as plt 
        [V_n, I_0] = Estimate_Node_voltage(A, Z, [i*1000 for i in opt.x], I0, V0, show_plot = True)       
                
        fig, axs = plt.subplots(2, sharex=True)
    #    # Current
    #    #axs[0].plot(I0, color = 'b', label = 'I0')
    #    #axs[0].plot(opt.x, color = 'r', label = 'Opt.x')
    #    #axs[0].set_ylabel('Node current, $I_{n}$, [A]')
    #    #axs[0].legend(loc='upper right')
    #    
    #    # Power    
        axs[0].plot(Sn0, color = 'b', label = 'Sn0')
        axs[0].plot([i for i in opt.x], color = 'r', label = 'Opt.x')
        axs[0].set_ylabel('Node power, $S_{n}$, [kW]')
        axs[0].legend(loc='upper right')
        axs[0].grid()    
    #    
    #    
    #    # Voltage
        axs[1].plot([i/400 for i in V0], color = 'b', label = 'V0')
        #axs[1].plot([i/400 for i in matmul(A,V0) - matmul(Z,(matmul(transpose(A),transpose(I0)) + divide(Sn0, V0))) + B*400], color = 'b', label = 'V0')
        #axs[1].plot([i/400 for i in matmul(A,V0) - matmul(Z,opt.x) + B*400], color = 'r', label = 'Opt.x')
        axs[1].plot([i/400 for i in V_n], color = 'r', label = 'Opt.x')
        axs[1].set_ylabel('Node voltage, $V_{n}$, [p.u.]')
        axs[1].legend(loc='upper right', ncol = 2)
        axs[1].set_xlim([0,301])
        axs[1].grid() 
        axs[1].set_xlabel('Node')     
         
    
    #    axs[2].grid()    
    #    axs[2].plot([P-S for S,P in zip(Sn0, [i for i in opt.x])])
    #    axs[2].set_ylabel('Change in node power, $\Delta S_{n}$, [kW]')


        
    return opt.x


#import pandas as pd
#from numpy import zeros
#DF_Network = pd.read_excel('Stedin_network_3.xlsx') # ['Gaia_network_2.xlsx', 'Stedin_network_2.xlsx', 'Stedin_network_3']
#Network_headers = ['Van.Nummer', 'Naar.Nummer', 'Lengte', 'Kabeltype', 'GM', 'House_Type']
#Wires_headers = ['Wire', 'R', 'X'] 
#wire_file = 'Gaia_cables.xlsx'
#csv_name = 'S_n_Stedin_301.csv'
#
#
## Define the arguments
#ts = 40
#
#Sn0 = Case_0_base_W[0][12].loc[1].tolist()
#V0 = Case_0_base_W[0][27].loc[1].tolist()
#I0 = Case_0_base_W[0][28].loc[1].tolist()
#
#Sn1 = Case_0_base_W[0][12].loc[ts].tolist()
#V1 = Case_0_base_W[0][27].loc[ts].tolist()
#I1 = Case_0_base_W[0][28].loc[ts].tolist()
#
#Sn4 = Case_4_MCES_NA_S[-1][12].loc[ts].tolist()
#V4 = Case_4_MCES_NA_S[-1][27].loc[ts].tolist()
#I4 = Case_4_MCES_NA_S[-1][28].loc[ts].tolist()
#
#Sn5 = Case_4_MCES_NA_S[-1][12].loc[340].tolist()
#V5 = Case_4_MCES_NA_S[-1][27].loc[340].tolist()
#I5 = Case_4_MCES_NA_S[-1][28].loc[340].tolist()
#
#[A, Z] = Create_network(DF_Network, Network_headers, Wires_headers, wire_file = 'Gaia_cables.xlsx')
#B = zeros(len(I0))
#B[0] = 1
#
#
#pu_lim = 0.95
#
#opt_results_0 = Find_power_setpoints(A, Z, B, V0, I0, Sn0, all_nodes, selected_nodes_list[-1][0], pu_lim)
#
#opt_results_1 = Find_power_setpoints(A, Z, B, V1, I1, Sn1, all_nodes, selected_nodes_list[-1][0], pu_lim)
#
#opt_results_4 = Find_power_setpoints(A, Z, B, V4, I4, Sn4, all_nodes, selected_nodes_list[-1][0], pu_lim)
#
#opt_results_5 = Find_power_setpoints(A, Z, B, V5, I5, Sn5, all_nodes, selected_nodes_list[-1][0], pu_lim)

###############################################################################
    
def Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True, prices = [1.15, 6500, 10000, 25000], module_power_ref = 400):
    from numpy import clip
    
    PV_sizes = Create_PV_DF(S_n, return_modules = True) 
    CAPEX = [0 for i in range(len(S_n[0]))]
    
    for node in range(len(S_n[0])):
        if node in selected_nodes:
            if PV:
                CAPEX[node] += clip(PV_sizes[node]*module_power_ref*prices[0], 2500, None)
            if HP:
                CAPEX[node] += prices[1]
            if BESS:
                CAPEX[node] += prices[2]
            if TESS:
                CAPEX[node] += prices[3]

    return [CAPEX[i] for i in selected_nodes]
            
###############################################################################
def Calculate_EoL(Case, C_BoL = 10):
    
    EoL = [0 for i in range(len(Case[0][4].tolist()))]
    for node in Case[0][3]:
#        print(node)
        EoL[node] = 0.2*C_BoL/Case[0][16][str(node+1)].tolist()[-1]
    
    return EoL

###############################################################################
    
def Calculate_Revenue(Case_n, Case_0, selected_nodes, case = 0, electricity_price = False, gas_price = 1.44, dt = 0.25):
    
#    Cost_base = [0 for i in range(len(Case_n[case][4]))]
#    Cost_n = [0 for i in range(len(Case_n[case][4]))]
    Cost_base = []
    Cost_n = []
    
    if not electricity_price:
        for node in range(len(Case_n[case][4])):
            if node in selected_nodes:
    #            Cost_base[node] = sum([p*c for p,c in zip(Case_0[case][12][str(node+1)].tolist(), Case_0[case][29]['Current'].tolist())])
    #            Cost_n[node] = sum([p*c for p,c in zip(Case_n[case][12][str(node+1)].tolist(), Case_n[case][29]['Current'].tolist())])
                Cost_base.append(sum([p*c for p,c in zip(Case_0[case][12][str(node+1)].tolist(), Case_0[case][29]['Current'].tolist())]) + gas_price*dt*Case_0[0][23][str(node+1)].sum()/(1000*12.7))
                Cost_n.append(sum([p*c for p,c in zip(Case_n[case][12][str(node+1)].tolist(), Case_n[case][29]['Current'].tolist())]) + gas_price*dt*Case_n[0][23][str(node+1)].sum()/(1000*12.7))

    else:
        for node in range(len(Case_n[case][4])):
            if node in selected_nodes:
    #            Cost_base[node] = sum([p*c for p,c in zip(Case_0[case][12][str(node+1)].tolist(), Case_0[case][29]['Current'].tolist())])
    #            Cost_n[node] = sum([p*c for p,c in zip(Case_n[case][12][str(node+1)].tolist(), Case_n[case][29]['Current'].tolist())])
                Cost_base.append(sum([p*electricity_price for p in Case_0[case][12][str(node+1)].tolist()]) + gas_price*dt*Case_0[0][23][str(node+1)].sum()/(1000*12.7))
                Cost_n.append(sum([p*electricity_price for p in Case_n[case][12][str(node+1)].tolist()]) + gas_price*dt*Case_n[0][23][str(node+1)].sum()/(1000*12.7))
        
        
#    return [Cost_base, Cost_n]
    
    return [base-n for base,n in zip(Cost_base, Cost_n)]

###############################################################################
    
#def Calculate_ROI(Case_n, Case_0, selected_nodes, PV = True, HP = True, BESS = True, TESS = True):
#    
#    Caluclate_CAPEX(Case_n[0][12].tolist(), selected_nodes, PV, HP, BESS, TESS)
#    
#    Calculate_Revenue(Case_1, Case_0, selected_nodes)

###############################################################################

#Plot_ROI(all_nodes, S_n, Case_0_base_Y, Case_1_PV_Y, Case_2_PV_HP_Y, Case_3_PV_HP_BESS_Y, Case_5 = Case_5_MCES_PA_Y)

def Plot_ROI(selected_nodes, S_n, Case_0, Case_1 = False, Case_2 = False, Case_3 = False, Case_4 = False, Case_5 = False, Case_6 = False, Case_7 = False, gas_price = 1.44, electricity_price = False, dt = 0.25):
    import matplotlib.pylab as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })

    
    E_load_ref = [(365/7)*dt*sum([i for i in Case_0[0][8][str(node+1)].tolist()[0:7*96]]) for node in selected_nodes]
    E_load_ref = [i/1000 for i in E_load_ref]


    plt.figure(constrained_layout=True)
    plt.grid()
    plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
    plt.xlim((0, max(E_load_ref)+1))
    plt.title('CAPEX for each case')
    plt.xlabel('Node yearly consumption [MWh]')
    plt.ylabel('CAPEX [€]')
    
    if Case_1:
        plt.scatter(E_load_ref, Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = False, BESS = False, TESS = False), label = 'Case 1')
    if Case_2:
        plt.scatter(E_load_ref, Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = False, TESS = False), color = 'r', label = 'Case 2')
    if Case_3:
        plt.scatter(E_load_ref, Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = False), color = 'b', label = 'Case 3')
    if Case_7:
        plt.scatter(E_load_ref, Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = False), color = 'b', label = 'Case 4')        
    if Case_4:
        plt.scatter(E_load_ref, Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True), color = 'y', label = 'Case 5')
    if Case_5:
        plt.scatter(E_load_ref, Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True), color = 'k', label = 'Case 6')
    if Case_6:
        plt.scatter(E_load_ref, Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True), color = 'g', label = 'Case 7')
    
    plt.legend(loc='lower right', ncol = 1)
    plt.show()     
    

    if not electricity_price:
        plt.figure(constrained_layout=True)
        plt.grid()
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('Revenue for each case')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('Revenue [€/year]')
        
        if Case_1:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_1, Case_0, selected_nodes), label = 'Case 1')
        if Case_2:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_2, Case_0, selected_nodes), color = 'r', label = 'Case 2')
        if Case_3:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_3, Case_0, selected_nodes), color = 'b', label = 'Case 3')
        if Case_7:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_7, Case_0, selected_nodes), color = 'b', label = 'Case 4')                 
        if Case_4:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_4, Case_0, selected_nodes), color = 'y', label = 'Case 5')
        if Case_5:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_5, Case_0, selected_nodes), color = 'k', label = 'Case 6')
        if Case_6:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_6, Case_0, selected_nodes), color = 'g', label = 'Case 7')    
            
        plt.legend(loc='lower right', ncol = 1)
        plt.show()  
    
    
        plt.figure(constrained_layout=True)
        plt.grid()
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('ROI for each case')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('ROI [years]')
        
        if Case_1:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = False, BESS = False, TESS = False), Calculate_Revenue(Case_1, Case_0, selected_nodes))], label = 'Case 1')
        if Case_2:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = False, TESS = False), Calculate_Revenue(Case_2, Case_0, selected_nodes))], color = 'r', label = 'Case 2')
        if Case_3:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = False), Calculate_Revenue(Case_3, Case_0, selected_nodes))], color = 'b', label = 'Case 3')
        if Case_7:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = False), Calculate_Revenue(Case_7, Case_0, selected_nodes))], color = 'b', label = 'Case 4')
        if Case_4:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True), Calculate_Revenue(Case_4, Case_0, selected_nodes))], color = 'y', label = 'Case 5')
        if Case_5:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True), Calculate_Revenue(Case_5, Case_0, selected_nodes))], color = 'k', label = 'Case 6')
        if Case_6:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True), Calculate_Revenue(Case_6, Case_0, selected_nodes))], color = 'g', label = 'Case 7')
        
        plt.legend(loc='upper right', ncol = 1)
        plt.show()     

    else:
        plt.figure(constrained_layout=True)
        plt.grid()
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('Revenue for each case')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('Revenue [€/year]')
        
        if Case_1:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_1, Case_0, selected_nodes, electricity_price = electricity_price), label = 'Case 1')
        if Case_2:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_2, Case_0, selected_nodes, electricity_price = electricity_price), color = 'r', label = 'Case 2')
        if Case_3:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_3, Case_0, selected_nodes, electricity_price = electricity_price), color = 'b', label = 'Case 3')
        if Case_7:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_7, Case_0, selected_nodes, electricity_price = electricity_price), color = 'b', label = 'Case 4')
        if Case_4:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_4, Case_0, selected_nodes, electricity_price = electricity_price), color = 'y', label = 'Case 5')
        if Case_5:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_5, Case_0, selected_nodes, electricity_price = electricity_price), color = 'k', label = 'Case 6')
        if Case_6:
            plt.scatter(E_load_ref, Calculate_Revenue(Case_6, Case_0, selected_nodes, electricity_price = electricity_price), color = 'g', label = 'Case 7')
          
        plt.legend(loc='lower right', ncol = 1)
        plt.show()  
    
    
        plt.figure(constrained_layout=True)
        plt.grid()
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('ROI for each case')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('ROI [years]')
        
        if Case_1:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = False, BESS = False, TESS = False), Calculate_Revenue(Case_1, Case_0, selected_nodes, electricity_price = electricity_price))], label = 'Case 1')
        if Case_2:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = False, TESS = False), Calculate_Revenue(Case_2, Case_0, selected_nodes, electricity_price = electricity_price))], color = 'r', label = 'Case 2')
        if Case_3:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = False), Calculate_Revenue(Case_3, Case_0, selected_nodes, electricity_price = electricity_price))], color = 'b', label = 'Case 3')
        if Case_7:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = False), Calculate_Revenue(Case_7, Case_0, selected_nodes, electricity_price = electricity_price))], color = 'b', label = 'Case 4')
        if Case_4:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True), Calculate_Revenue(Case_4, Case_0, selected_nodes, electricity_price = electricity_price))], color = 'y', label = 'Case 5')
        if Case_5:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True), Calculate_Revenue(Case_5, Case_0, selected_nodes, electricity_price = electricity_price))], color = 'k', label = 'Case 6')
        if Case_6:
            plt.scatter(E_load_ref, [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True), Calculate_Revenue(Case_6, Case_0, selected_nodes, electricity_price = electricity_price))], color = 'g', label = 'Case 7')
        
        plt.legend(loc='upper right', ncol = 1)
        plt.show()


    if Case_5:
        ROI_2 = [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = False, TESS = False), Calculate_Revenue(Case_2, Case_0, selected_nodes))]
        ROI_7 = [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = False), Calculate_Revenue(Case_7, Case_0, selected_nodes))]
        ROI_5 = [C/R for C,R in zip(Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True), Calculate_Revenue(Case_5, Case_0, selected_nodes))]
        
        
        CAPEX_2 = Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = False, TESS = False)
        CAPEX_7 = Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = False)
        CAPEX_5 = Caluclate_CAPEX(S_n, selected_nodes, PV = True, HP = True, BESS = True, TESS = True)
        
        
        Revenue_2 = Calculate_Revenue(Case_2, Case_0, selected_nodes)
        Revenue_7 = Calculate_Revenue(Case_7, Case_0, selected_nodes)
        Revenue_5 = Calculate_Revenue(Case_5, Case_0, selected_nodes)
        
        Grid_Cost_2 = [sum([p*c for p,c in zip(Case_2[0][12][str(node+1)].tolist(), Case_2[0][29]['Current'].tolist())]) + gas_price*dt*Case_2[0][23][str(node+1)].sum()/(1000*12.7) for node in selected_nodes]
        Grid_Cost_7 = [sum([p*c for p,c in zip(Case_7[0][12][str(node+1)].tolist(), Case_2[0][29]['Current'].tolist())]) + gas_price*dt*Case_7[0][23][str(node+1)].sum()/(1000*12.7) for node in selected_nodes]
        Grid_Cost_5 = [sum([p*c for p,c in zip(Case_5[0][12][str(node+1)].tolist(), Case_5[0][29]['Current'].tolist())]) + gas_price*dt*Case_5[0][23][str(node+1)].sum()/(1000*12.7) for node in selected_nodes]    

        Compensation_7 = [Calculate_Compensation(c2,c5,r2,r5) for c2,c5,r2,r5 in zip(CAPEX_2, CAPEX_7, Revenue_2, Revenue_7)]
        Compensation_5 = [Calculate_Compensation(c2,c5,r2,r5) for c2,c5,r2,r5 in zip(CAPEX_2, CAPEX_5, Revenue_2, Revenue_5)]
#        
#        print('ROI_2: ', ROI_2[0])
#        print('ROI_5: ', ROI_5[0])
#        
#        print('CAPEX_2: ', CAPEX_2[0])
#        print('CAPEX_5: ', CAPEX_5[0])
#        
#        print('Revenue_2: ', Revenue_2[0])
#        print('Revenue_5: ', Revenue_5[0])        
#        

        plt.figure(constrained_layout=True)
        plt.grid()
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('ROI ratio for Cases 2 and 4')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('ROI ratio [-]')
        plt.scatter(E_load_ref, [R7/R2 for R2,R7 in zip(ROI_2, ROI_7)])
        plt.show()


        plt.figure(constrained_layout=True)
        plt.grid()
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('Compensation required for case 4')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('[€/month]')
        plt.scatter(E_load_ref, Compensation_7, label='Compensation required')
        plt.scatter(E_load_ref, [i/12 for i in Grid_Cost_7], color = 'r', label='Average energy cost')
        plt.legend(loc='lower right')
        plt.show()        

        plt.figure(constrained_layout=True)
        plt.grid()
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('Difference in compensation and cost for case 4')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('Compensation - Cost, [€/month]')
        plt.scatter(E_load_ref, [cost - compensation for cost,compensation in zip ([i/12 for i in Grid_Cost_7], Compensation_7)], label='Compensation required')
        plt.show() 

        plt.figure(constrained_layout=True)
        plt.grid()
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('ROI ratio for Cases 2 and 6')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('ROI ratio [-]')
        plt.scatter(E_load_ref, [R5/R2 for R2,R5 in zip(ROI_2, ROI_5)])
        plt.show()


        plt.figure(constrained_layout=True)
        plt.grid()
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('Compensation required for case 6')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('[€/month]')
        plt.scatter(E_load_ref, Compensation_5, label='Compensation required')
        plt.scatter(E_load_ref, [i/12 for i in Grid_Cost_5], color = 'r', label='Average energy cost')
        plt.legend(loc='lower right')
        plt.show()        

        plt.figure(constrained_layout=True)
        plt.grid()
        plt.xticks(arange(0, max(E_load_ref)+2, step = 1))
        plt.xlim((0, max(E_load_ref)+1))
        plt.title('Difference in compensation and cost for case 6')
        plt.xlabel('Node yearly consumption [MWh]')
        plt.ylabel('Compensation - Cost, [€/month]')
        plt.scatter(E_load_ref, [cost - compensation for cost,compensation in zip ([i/12 for i in Grid_Cost_5], Compensation_5)], label='Compensation required')
        plt.show()  

###############################################################################
        
def Calculate_Compensation(CAPEX_2, CAPEX_5, Revenue_2, Revenue_5):
    
    return ((CAPEX_5/CAPEX_2)*Revenue_2 - Revenue_5)/12

###############################################################################
    
def define_worse_nodes(case_n, n_nodes = 301, n_cases = 5, n_timesteps = 7*96, v_feeder = 400):
    
    n_violations_20 = [0 for i in range(n_nodes)]
    n_violations_40 = [0 for i in range(n_nodes)]
    n_violations_60 = [0 for i in range(n_nodes)]
    n_violations_80 = [0 for i in range(n_nodes)]
    
    sum_voltage_incompliance_20 = [0 for i in range(n_nodes)]
    sum_voltage_incompliance_40 = [0 for i in range(n_nodes)]
    sum_voltage_incompliance_60 = [0 for i in range(n_nodes)]
    sum_voltage_incompliance_80 = [0 for i in range(n_nodes)]
    
    
    for node in range(n_nodes):
        for case in range(n_cases):
            for ts in range(n_timesteps):
                if case_n[0+case][4][node][ts]/v_feeder < 0.95:
                    n_violations_20[node] += 1
                    sum_voltage_incompliance_20[node] += 0.95 - case_n[0+case][4][node][ts]/v_feeder
                elif case_n[0+case][4][node][ts]/v_feeder > 1.05:
                    n_violations_20[node] += 1
                    sum_voltage_incompliance_20[node] += case_n[0+case][4][node][ts]/v_feeder - 1.05

                if case_n[5+case][4][node][ts]/v_feeder < 0.95:
                    n_violations_40[node] += 1
                    sum_voltage_incompliance_40[node] += 0.95 - case_n[5+case][4][node][ts]/v_feeder
                elif case_n[5+case][4][node][ts]/v_feeder > 1.05:
                    n_violations_40[node] += 1
                    sum_voltage_incompliance_40[node] += case_n[5+case][4][node][ts]/v_feeder - 1.05

                if case_n[10+case][4][node][ts]/v_feeder < 0.95:
                    n_violations_60[node] += 1
                    sum_voltage_incompliance_60[node] += 0.95 - case_n[10+case][4][node][ts]/v_feeder
                elif case_n[10+case][4][node][ts]/v_feeder > 1.05:
                    n_violations_60[node] += 1
                    sum_voltage_incompliance_60[node] += case_n[10+case][4][node][ts]/v_feeder - 1.05
                    
                if case_n[15+case][4][node][ts]/v_feeder < 0.95:
                    n_violations_80[node] += 1
                    sum_voltage_incompliance_80[node] += 0.95 - case_n[15+case][4][node][ts]/v_feeder
                elif case_n[15+case][4][node][ts]/v_feeder > 1.05:
                    n_violations_80[node] += 1
                    sum_voltage_incompliance_80[node] += case_n[15+case][4][node][ts]/v_feeder - 1.05                

    import matplotlib.pylab as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })
    
    
#    plt.figure(constrained_layout=True)
#    plt.grid()
#    plt.plot([i for i in range(n_nodes)], n_violations_20, label = '20 %')
#    plt.plot([i for i in range(n_nodes)], n_violations_40, color = 'r', label = '40 %')
#    plt.plot([i for i in range(n_nodes)], n_violations_60, color = 'b', label = '60 %')
#    plt.plot([i for i in range(n_nodes)], n_violations_80, color = 'y', label = '80 %')
#    plt.xlim([0,n_nodes])
#    plt.xlabel('Node')
#    plt.ylabel('Number of voltage incompliance')
#    plt.legend()
#    plt.show()
#
#
#    plt.figure(constrained_layout=True)
#    plt.grid()
#    plt.plot([i for i in range(n_nodes)], sum_voltage_incompliance_20, label = '20 %')
#    plt.plot([i for i in range(n_nodes)], sum_voltage_incompliance_40, color = 'r', label = '40 %')
#    plt.plot([i for i in range(n_nodes)], sum_voltage_incompliance_60, color = 'b', label = '60 %')
#    plt.plot([i for i in range(n_nodes)], sum_voltage_incompliance_80, color = 'y', label = '80 %')
#    plt.xlim([0,n_nodes])
#    plt.xlabel('Node')
#    plt.ylabel('Sum of voltage incompliance')
#    plt.legend()
#    plt.show()    


    from numpy import array
    [n_violations_20, n_violations_40, n_violations_60, n_violations_80] =  [array(n_violations_20), array(n_violations_40), array(n_violations_60), array(n_violations_80)]
    [sum_voltage_incompliance_20, sum_voltage_incompliance_40, sum_voltage_incompliance_60, sum_voltage_incompliance_80] = [array(sum_voltage_incompliance_20), array(sum_voltage_incompliance_40), array(sum_voltage_incompliance_60), array(sum_voltage_incompliance_80)]
    

#    plt.figure(constrained_layout=True)
#    plt.grid()
#    plt.bar(array([i for i in range(n_nodes)]), n_violations_20, label = '20 %')
#    plt.bar(array([i for i in range(n_nodes)]), n_violations_40, color = 'r', label = '40 %', bottom = n_violations_20)
#    plt.bar(array([i for i in range(n_nodes)]), n_violations_60, color = 'b', label = '60 %', bottom = n_violations_20+n_violations_40)
#    plt.bar(array([i for i in range(n_nodes)]), n_violations_80, color = 'y', label = '80 %', bottom = n_violations_20+n_violations_40+n_violations_60)
#    plt.xlim([0, n_nodes])
#    plt.xlabel('Node')
#    plt.ylabel('Number of voltage incompliance')
#    plt.legend(ncol=2, prop={'size': 12})
#    plt.show()


    plt.figure(constrained_layout=True)
    plt.grid()
    plt.bar(array([i for i in range(n_nodes)]), sum_voltage_incompliance_20, label = '20 %')
    plt.bar(array([i for i in range(n_nodes)]), sum_voltage_incompliance_40, color = 'r', label = '40 %', bottom = sum_voltage_incompliance_20)
    plt.bar(array([i for i in range(n_nodes)]), sum_voltage_incompliance_60, color = 'b', label = '60 %', bottom = sum_voltage_incompliance_20+sum_voltage_incompliance_40)
    plt.bar(array([i for i in range(n_nodes)]), sum_voltage_incompliance_80, color = 'y', label = '80 %', bottom = sum_voltage_incompliance_20+sum_voltage_incompliance_40+sum_voltage_incompliance_60)
    plt.xlim([0, n_nodes])
    plt.xlabel('Node')
    plt.ylabel('Sum of voltage incompliance. $\Sigma_{V, i}$')
    plt.ylim([0, 100])
    plt.legend(ncol=2, prop={'size': 12})
    plt.show()   
    
    
###############################################################################
    
def Calculate_P_BESS_compensation(case_n, case, ts, n_nodes = 301):
    
    return case_n[case][12].iloc[ts].sum()

###############################################################################
    
def Size_BESS(case_n, case, n_timesteps = 7*96, start_day = 0, P_PV = 330, v_feeder = 400, module_power_ref = 0.315, dt = 0.25):     # day_ranges = 177
    
    P_BESS = [0 for i in range(n_timesteps)]
    V_min = []
    
    for ts in range(n_timesteps):
#        print(ts)
        V_min.append(min([i[ts] for i in case_n[0][4]]))
        if (min([i[ts] for i in case_n[0][4]]) < 0.95*v_feeder) or (max([i[ts] for i in case_n[0][4]]) > 1.05*v_feeder):
            P_BESS[ts] = Calculate_P_BESS_compensation(case_n, case, ts)

    
    n_days = int(n_timesteps/(24/dt))        
    daily_energy_BESS = [dt*sum(P_BESS[int(day*24/dt):int((day+1)*24/dt)]) for day in range(n_days)]
    
    P_charge_BESS = [daily_energy_BESS[day]/(dt*P_BESS[int(day*24/dt):int((day+1)*24/dt)].count(0)) for day in range(n_days)]

#    print(P_charge_BESS)

    
    
    P_BESS = [P if P != 0 else -P_charge_BESS[int(ts//(24/dt))] for P,ts in zip(P_BESS,range(int(n_days*24/dt)))]

    import csvreader
    
    CSVDataPV = csvreader.read_data(csv='PV_15min.csv', address='')
    CSVDataPV.data2array()

    Reference_PV_profile = [(1/module_power_ref)*i[0]/1000 for i in CSVDataPV.ar.tolist()]    # 1kWp
    daily_energy_PV = [dt*sum(Reference_PV_profile[int((start_day + day)*24/dt):int((start_day + day+1)*24/dt)]) for day in range(n_days)]

    PV_sizes = [BESS/PV for BESS,PV in zip(daily_energy_BESS, daily_energy_PV)]
        
        

    import matplotlib.pylab as plt

    plt.rcParams.update({
    #    "text.usetex": True,
        "font.family": "Times New Roman",
        'font.size': 16
    })

    plt.figure(constrained_layout=True)
    plt.grid()
    plt.plot(V_min)
    plt.xlim([0, n_timesteps])
    plt.xlabel('Time')
    plt.ylabel('Minimum node voltage, $V_{min}$, [V]')
    plt.title('Minimum node voltage')
        
    
    plt.figure(constrained_layout=True)
    plt.grid()
    plt.plot(P_BESS)
    plt.xlim([0, n_timesteps])
    plt.xlabel('Time')
    plt.ylabel('BESS power, $P_{BESS}$, [kW]')
    plt.title('Centralized BESS power')

    plt.figure(constrained_layout=True)
    plt.grid()
    plt.bar([i+1 for i in range(n_days)], daily_energy_BESS)
    plt.xlim([0, n_days+1])
    plt.xlabel('Time')
    plt.ylabel('Daily BESS energy, $E_{BESS}/day$, [kWh]')
    plt.title('Daily BESS energy')

    plt.figure(constrained_layout=True)
    plt.grid()
    plt.bar([i+1 for i in range(n_days)], daily_energy_PV)
    plt.xlim([0, n_days+1])
    plt.xlabel('Time')
    plt.ylabel('Daily PV energy, $E_{PV}/day$, [kWh]')
    plt.title('Daily PV energy')    

    plt.figure(constrained_layout=True)
    plt.grid()
    plt.bar([i+1 for i in range(n_days)], PV_sizes)
    plt.xlim([0, n_days+1])
    plt.xlabel('Time')
    plt.ylabel('Daily PV size, $P_{PV}$, [kWp]')
    plt.title('Daily PV size')    
    
    return P_BESS


###############################################################################
################################   Main code   ################################
    

#%%############################################################################
################################   CIGRE-18   #################################
global t_registry
t_registry = []
  
DF_Network = pd.read_excel('Cigre_18_2.xlsx') 
Network_headers = ['Van.Nummer', 'Naar.Nummer', 'Lengte', 'Kabeltype', 'GM', 'House_Type'] # ['From', 'To', 'Distance', 'Wire', 'GM', 'House_Type']

wire_file = 'Wires.xlsx'
Wires_headers = ['Wires', 'resistance', 'impedance'] 

csv_name = 'S_n_CIGRE_18.csv' 

season = 'Test'
penetration = 100
case = 0

S_n = 3*Create_loads(DF_Network, Network_headers, Load_data = True, csv_name = csv_name)/1000
#S_n = Create_loads(DF_Network, Network_headers)

selected_nodes = [10, 14, 15, 16, 17]
remaining_nodes = []
all_nodes = [i for i in range(18)]

[enable_HP, BESS_power, BESS_energy, enable_TESS, TESS_capacity, follow_control] = case_parameters(3)
force_temperature = True


control_type = 'heuristic'


if season == 'Winter':
    start_day = 0# [0, 177]       # Winter = 0, Summer = 177,
    end_day =  7#[7, 184]         # Winter = 7, Summer = 184,

elif season == 'Summer':
    start_day = 177# [0, 177]       # Winter = 0, Summer = 177,
    end_day =  184#[7, 184]         # Winter = 7, Summer = 184,

elif season == 'Year':
    start_day = 0# [0, 177]       # Winter = 0, Summer = 177,
    end_day =  364#[7, 184]         # Winter = 7, Summer = 184,  
    
elif season == 'Test':
    start_day = 0# [0, 177]       # Winter = 0, Summer = 177,
    end_day =  1#[7, 184]         # Winter = 7, Summer = 184,    
    

dt = 0.25
t_final = int((end_day - start_day)*24/dt)    

Case_CIGRE = []
#Case_CIGRE_controlled = []

#A = Create_admittance_matrix(DF_Network, Network_headers, show_heatmap = True)
#Z = Create_impedance_matrix(DF_Network, Network_headers, Wires_headers, show_heatmap = True, wire_file = wire_file, km = True)


start = time.time()    # The timer is initializad.

[V_n_registry, I_registry, [df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in, df_P_n, df_V_n, df_I_n, df_Prices, df_House_type]] = simulate_Network(DF_Network, Network_headers, Wires_headers, selected_nodes, remaining_nodes, S_n, t_final, t_0 = int(24*start_day/dt), control_type = control_type, wire_file = wire_file, follow_control = False, P_BESS_max = BESS_power, Capacity_BESS_BoL = BESS_energy, enable_HP = enable_HP, enable_TESS = enable_TESS, force_temperature = force_temperature)
Case_CIGRE.append([season, penetration, case, selected_nodes, V_n_registry, I_registry, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in, df_P_n, df_V_n, df_I_n, df_Prices, df_House_type])

#[V_n_registry_controlled, I_registry_controlled, [df_P_PV_av_controlled, df_P_PV_controlled, df_P_L_controlled, df_P_HP_controlled, df_P_HP_TESS_controlled, df_P_BESS_controlled, df_P_Grid_controlled, df_P_ref_controlled, df_DSO_operation_controlled, df_SOC_BESS_controlled, df_C_BESS_controlled, df_Qdot_D_controlled, df_Qdot_HP_controlled, df_Qdot_HP_TESS_controlled, df_Qdot_TESS_controlled, df_T_TESS_controlled, df_Qdot_TESS_SD_controlled, df_Qdot_Boiler_controlled, df_enable_HP_2_TESS_controlled, df_T_in_controlled, df_P_n_controlled, df_V_n_controlled, df_I_n_controlled, df_Prices_controlled, df_House_type_controlled]] = simulate_Network(DF_Network, Network_headers, Wires_headers, selected_nodes, remaining_nodes, S_n, t_final, t_0 = int(24*start_day/dt), control_type = control_type, wire_file = wire_file, follow_control = follow_control, P_BESS_max = BESS_power, Capacity_BESS_BoL = BESS_energy, enable_HP = enable_HP, enable_TESS = enable_TESS)
#Case_CIGRE_controlled.append([season, penetration, case, selected_nodes, V_n_registry_controlled, I_registry_controlled, df_P_PV_av_controlled, df_P_PV_controlled, df_P_L_controlled, df_P_HP_controlled, df_P_HP_TESS_controlled, df_P_BESS_controlled, df_P_Grid_controlled, df_P_ref_controlled, df_DSO_operation_controlled, df_SOC_BESS_controlled, df_C_BESS_controlled, df_Qdot_D_controlled, df_Qdot_HP_controlled, df_Qdot_HP_TESS_controlled, df_Qdot_TESS_controlled, df_T_TESS_controlled, df_Qdot_TESS_SD_controlled, df_Qdot_Boiler_controlled, df_enable_HP_2_TESS_controlled, df_T_in_controlled, df_P_n_controlled, df_V_n_controlled, df_I_n_controlled, df_Prices_controlled, df_House_type_controlled])

end = time.time()    # The timer is initializad.
totalelapsed = end - start  # The total time is calculated.


Make_voltage_boxplots([V_n_registry], ['CIGRE 18'])
#Voltage_Power_correlation(all_nodes, df_P_Grid, df_V_n, df_Prices, tf = t_final)

for node in selected_nodes:

    Plot_house_behaviour(Case_CIGRE, node+1, end = t_final) # Apartment