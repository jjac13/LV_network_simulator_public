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



#%%############################################################################
################################   CIGRE-18   #################################

import pandas as pd
from numpy import array
from Network_builder_public import *

global t_registry
t_registry = []
  
DF_Network = pd.read_excel('Cigre_18_2.xlsx') 
Network_headers = ['Van.Nummer', 'Naar.Nummer', 'Lengte', 'Kabeltype', 'GM', 'House_Type']

wire_file = 'Wires.xlsx'
Wires_headers = ['Wires', 'resistance', 'impedance'] 
csv_name = 'S_n_CIGRE_18.csv' 

season = 'Winter'


S_n = 3*Create_loads(DF_Network, Network_headers, Load_data = True, csv_name = csv_name)/1000


all_nodes = [i for i in range(18)]

force_temperature = True
control_type = 'heuristic'


if season == 'Winter':
    start_day = 0       # Winter = 0, Summer = 177,
    end_day =  7        # Winter = 7, Summer = 184,

elif season == 'Summer':
    start_day = 177     # Winter = 0, Summer = 177,
    end_day =  184      # Winter = 7, Summer = 184,

elif season == 'Year':
    start_day = 0       # Winter = 0, Summer = 177,
    end_day =  364      # Winter = 7, Summer = 184,  
    
elif season == 'Test':
    start_day = 0       # Winter = 0, Summer = 177,
    end_day =  1        # Winter = 7, Summer = 184,    
    

dt = 0.25
t_final = int((end_day - start_day)*24/dt)    

Case_CIGRE = []
Case_CIGRE_controlled = []
Case_CIGRE_centralized = []



print('Simulating base case (PV and heat pumps)')
case = 2
penetration = 100
selected_nodes = [10, 14, 15, 16, 17]
remaining_nodes = []
[enable_HP, BESS_power, BESS_energy, enable_TESS, TESS_capacity, follow_control] = case_parameters(case)
[V_n_registry, I_registry, [df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in, df_P_n, df_V_n, df_I_n, df_Prices, df_House_type]] = simulate_Network(DF_Network, Network_headers, Wires_headers, selected_nodes, remaining_nodes, S_n, t_final, t_0 = int(24*start_day/dt), control_type = control_type, wire_file = wire_file, follow_control = False, P_BESS_max = BESS_power, Capacity_BESS_BoL = BESS_energy, enable_HP = enable_HP, enable_TESS = enable_TESS, force_temperature = force_temperature, start_day = start_day, end_day = end_day)
Case_CIGRE.append([season, penetration, case, selected_nodes, V_n_registry, I_registry, df_P_PV_av, df_P_PV, df_P_L, df_P_HP, df_P_HP_TESS, df_P_BESS, df_P_Grid, df_P_ref, df_DSO_operation, df_SOC_BESS, df_C_BESS, df_Qdot_D, df_Qdot_HP, df_Qdot_HP_TESS, df_Qdot_TESS, df_T_TESS, df_Qdot_TESS_SD, df_Qdot_Boiler, df_enable_HP_2_TESS, df_T_in, df_P_n, df_V_n, df_I_n, df_Prices, df_House_type])



print('Simulating aggregation scenario')
case = 5
penetration = 100
selected_nodes = [10, 14, 15, 16, 17]
remaining_nodes = []
[enable_HP, BESS_power, BESS_energy, enable_TESS, TESS_capacity, follow_control] = case_parameters(case)
[V_n_registry_controlled, I_registry_controlled, [df_P_PV_av_controlled, df_P_PV_controlled, df_P_L_controlled, df_P_HP_controlled, df_P_HP_TESS_controlled, df_P_BESS_controlled, df_P_Grid_controlled, df_P_ref_controlled, df_DSO_operation_controlled, df_SOC_BESS_controlled, df_C_BESS_controlled, df_Qdot_D_controlled, df_Qdot_HP_controlled, df_Qdot_HP_TESS_controlled, df_Qdot_TESS_controlled, df_T_TESS_controlled, df_Qdot_TESS_SD_controlled, df_Qdot_Boiler_controlled, df_enable_HP_2_TESS_controlled, df_T_in_controlled, df_P_n_controlled, df_V_n_controlled, df_I_n_controlled, df_Prices_controlled, df_House_type_controlled]] = simulate_Network(DF_Network, Network_headers, Wires_headers, selected_nodes, remaining_nodes, S_n, t_final, t_0 = int(24*start_day/dt), control_type = control_type, wire_file = wire_file, follow_control = follow_control, P_BESS_max = BESS_power, Capacity_BESS_BoL = BESS_energy, enable_HP = enable_HP, enable_TESS = enable_TESS, start_day = start_day, end_day = end_day)
Case_CIGRE_controlled.append([season, penetration, case, selected_nodes, V_n_registry_controlled, I_registry_controlled, df_P_PV_av_controlled, df_P_PV_controlled, df_P_L_controlled, df_P_HP_controlled, df_P_HP_TESS_controlled, df_P_BESS_controlled, df_P_Grid_controlled, df_P_ref_controlled, df_DSO_operation_controlled, df_SOC_BESS_controlled, df_C_BESS_controlled, df_Qdot_D_controlled, df_Qdot_HP_controlled, df_Qdot_HP_TESS_controlled, df_Qdot_TESS_controlled, df_T_TESS_controlled, df_Qdot_TESS_SD_controlled, df_Qdot_Boiler_controlled, df_enable_HP_2_TESS_controlled, df_T_in_controlled, df_P_n_controlled, df_V_n_controlled, df_I_n_controlled, df_Prices_controlled, df_House_type_controlled])


print('Simulating centralized BESS scenario')
selected_nodes = [17]  
remaining_nodes = []
Capacity_BESS_BoL = 25
P_BESS_max = 20
P_PV_peak = 20
P_PV_peak = array([P_PV_peak if node in selected_nodes else 0 for node in range(len(S_n[0]))])

[V_n_registry_centralized, I_registry_centralized, [df_P_L_centralized, df_P_BESS_centralized, df_P_Grid_centralized, df_P_ref_centralized, df_SOC_BESS_centralized, df_C_BESS_centralized, df_P_PV_av_centralized, df_P_PV_centralized, df_P_n_centralized, df_V_n_centralized, df_I_n_centralized]] = simulate_Network_Centralized_BESS(DF_Network, Network_headers, Wires_headers, selected_nodes, remaining_nodes, S_n, P_PV_peak, df_P_L, t_final, t_0 = int(start_day*24/dt), P_BESS_max = P_BESS_max, Capacity_BESS_BoL = Capacity_BESS_BoL, wire_file = wire_file, V_ref_min = 0.98*400, V_ref_max = 1.02*400)
Case_CIGRE_centralized.append([season, penetration, case, selected_nodes, V_n_registry_centralized, I_registry_centralized, df_P_L_centralized, df_P_BESS_centralized, df_P_Grid_centralized, df_P_ref_centralized, df_SOC_BESS_centralized, df_C_BESS_centralized, df_P_PV_av_centralized, df_P_PV_centralized, df_P_n_centralized, df_V_n_centralized, df_I_n_centralized])



# Data visualization

Make_voltage_boxplots([V_n_registry, V_n_registry_controlled, V_n_registry_centralized], ['Base', 'Aggregation', 'Centralized'])

print('Voltage compliance for different scenarios [[above limit], [below limit]]')
print('Base (+- 0.01 pu): ', count_V_compliance(Case_CIGRE, lim = 0.02))
print('Aggregation (+- 0.01 pu): ', count_V_compliance(Case_CIGRE_controlled, lim = 0.02))
print('Centralized (+- 0.01 pu): ', count_V_compliance(Case_CIGRE_centralized, lim = 0.02))


Voltage_Power_correlation(all_nodes, df_P_Grid, df_V_n, df_Prices, tf = t_final)

for node in [10]: # selected_nodes

    Plot_house_behaviour(Case_CIGRE, node+1, end = t_final) # Apartment
    Plot_house_behaviour(Case_CIGRE_controlled, node+1, end = t_final) # Apartment  

Plot_grid_behaviour(Case_CIGRE, Case_CIGRE_centralized, selected_nodes[0], start = start_day, end = end_day)