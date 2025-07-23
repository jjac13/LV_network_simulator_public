# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a library with a class that defines devices that combined
## can create a synthetic load profile for electric power consumed by
## households, based on their probability of being activated at a certain
## momemnt.
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 1.0



###############################################################################
###########################   Imported modules   ##############################

from numpy import array

############################   Classes definition  ############################

## The class HouseholdDevice allows to create a device. with its maximum power,
## hourly probability and load profile during the timestep. The class also
## contains several appliances as examples, with typical values.

## For the class HouseholdDevice, the description of the parameters is the 
## next:
##  Power: holds the maximum power for the device in kW.
##  HourlyProbability: (array) holds the probability for the device to be
##      activated at a certain hour.
##  TimestepProfile: (array) holds the load profile for the device, during the
##      chosen time step.


class HouseholdDevice:
    def __init__(self, Power = 0, HourlyProbability = [], TimestepProfile = []):
        self.Power = Power                          # Maximum power of the device, in kW.
        self.HourlyProbability = HourlyProbability  # Probabilistic threshold for the device to start at a certain time.
        self.TimestepProfile = TimestepProfile      # Behaviour of the device, during the timestep, given every minute.
        
    def Fridge():       
        Fridge = HouseholdDevice()
        Fridge.Power = 0.080
        Fridge.HourlyProbability = array([0,0,0,0,0,0,0,0,0,0,0,0,100,100,0,0,0,0,0,0,0,0,0,0,0,0,100,100,0,0,0,0,0,100,100,0,0,0,0,0,100,0,0,0,0,0,0,0,0,100,100,0,0,0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,0,0,0,0,100,100,0,0,0,0,0,0,0,100,100,0,0,0,0,0,0,0,0,0,0,0])
        Fridge.TimestepProfile = array([Fridge.Power]*15)
        
        return Fridge

    def Router():
        Router = HouseholdDevice()
        Router.Power = 0.018
        Router.HourlyProbability = array([100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100])
        Router.TimestepProfile = array([Router.Power]*15)
        
        return Router
    
    def Dishwasher():        
        Dishwasher = HouseholdDevice()
        Dishwasher.Power = 1.2
        Dishwasher.HourlyProbability = array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,100,0,100,0,0,0,0,0,0,0,0])
        Dishwasher.TimestepProfile = array([Dishwasher.Power]*15)
        
        return Dishwasher
        
    def Lights():
        Lights = HouseholdDevice()
        Lights.Power = 0.030
        Lights.HourlyProbability = array([20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,40,40,40,40,40,40,90,90,90,90,90,90,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,30,60,60,60,60,60,60,95,95,95,95,95,95,95,95,95,95,95,95,95,40,40,40,40,40,40,40,40,20,20,20])
        Lights.TimestepProfile = array([Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100, Lights.Power*SetProbability()/100])
        
        return Lights
    
    def TV():
        TV = HouseholdDevice()
        TV.Power = 0.060
        TV.HourlyProbability = array([15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,40,40,40,40,40,40,25,25,25,25,25,25,25,25,25,25,40,40,50,50,50,50,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,40,40,40,40,40,40,40,40,40,60,60,60,60,60,60,60,60,60,60,60,25,25,25,25,25,25,25])
        TV.TimestepProfile = array([TV.Power]*15)
        
        return TV
        
    def Cell_charger():
        Cell_charger = HouseholdDevice()
        Cell_charger.Power = 0.015*2
        Cell_charger.HourlyProbability = array([40,40,40,40,40,40,40,40,40,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,50,50,50,50,50,50,80,80,80,80,80,80,50,50])
        Cell_charger.TimestepProfile = array([Cell_charger.Power]*15)
        
        return Cell_charger     
    
    def Laptop():
        Laptop = HouseholdDevice()
        Laptop.Power = 0.040*2
        Laptop.HourlyProbability = array([30,30,30,30,30,30,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,50,50,50,50,90,90,90,90,90,90,90,90,90,90,90,75,40,40,40,40,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,50,50,50,50,20,20,20,20,20,20,20,20,20,20,20,20,30,30,30,30,30,30,30,30,30,30,30,30])
        Laptop.TimestepProfile = array([Laptop.Power]*15)
        
        return Laptop
        
    def Screen():
        Screen = HouseholdDevice()
        Screen.Power = 0.040*2
        Screen.HourlyProbability = array([30,30,30,30,30,30,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,50,50,50,50,90,90,90,90,90,90,90,90,90,90,90,75,40,40,40,40,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,50,50,50,50,20,20,20,20,20,20,20,20,20,20,20,20,30,30,30,30,30,30,30,30,30,30,30,30])
        Screen.TimestepProfile = array([Screen.Power]*15)
        
        return Screen
        
    def Kitchen():
        Kitchen = HouseholdDevice()
        Kitchen.Power = 1.2
        Kitchen.HourlyProbability = array([15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,60,60,60,60,40,40,25,25,25,25,25,25,25,25,25,25,25,25,60,60,60,50,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,80,80,50,50,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20])
        Kitchen.TimestepProfile = array([Kitchen.Power]*10 + [0]*5)
        
        return Kitchen
    
    def Microwave():        
        Microwave = HouseholdDevice()
        Microwave.Power = 1.2
        Microwave.HourlyProbability = array([30,30,30,30,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,50,50,50,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,50,50,15,15,15,15,15,15,15,15,15,15,15,15,15,30,30,15,15,15,50,15,15,30,30])
        Microwave.TimestepProfile = array([Microwave.Power]*2 + [0]*13)
        
        return Microwave
        
    def Kettle():
        Kettle = HouseholdDevice()
        Kettle.Power = 1.5
        Kettle.HourlyProbability = array([15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,40,90,90,60,60,40,15,15,15,15,50,50,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,40,40,60,60,60,60,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,50,50,15,15,15,15,15,15,15,15,15,15,15,15])
        Kettle.TimestepProfile = array([Kettle.Power]*3 + [0]*12)    
        
        return Kettle
    
    
## The function SetProbability returns a random number between 0 and 100. The
## probability of each value has a normal distribution, with the higher numbers
## having higher probabilities of be returned.
        
def SetProbability():
    from numpy import random
    Probability = abs(random.standard_normal())
    if Probability >= 3:
        while Probability >= 3:
            Probability = abs(random.standard_normal())
    return 100-Probability*(100/3)


## The function SetBackgroundLoad creates a background load profile, based on 
## the devices indicated, during the time elapsed between the start and end
## times, with the timestep indicated.

## For the function SetBackgroundLoad, the description of the parameters is the 
## next:
##  BackgroundDevices: (list) contains an array with the devices considered
##      for the background load. Each element should be a HouseholdDevice.
##  EndTime: holds the end time for the simulation, in hours.
##  TimeStep: holds the sampling time for the resulting background load profile
##      in hours. It is considered as default a 15 min (0.25 h) step.
##  StartTime: holds the end time for the simulation, in hours. It is 
##      considered that the simulation starts at 0 time.
        
## Considerations:
##  - The resulting array will have as minimum timestep the timestep used for
##    the device load profile.

def SetBackgroundLoad(BackgroundDevices, EndTime, TimeStep = 0.25, StartTime = 0):
    from numpy import append
    
    BackgroundLoad = array([0]*60*int((EndTime-StartTime)))

    for device in BackgroundDevices:
        devicePower = array([])       
    
        for instant in range(int(len(BackgroundLoad)/TimeStep/60)):
            devicePower = append(devicePower, array([i*device.HourlyProbability[instant]/100 for i in device.TimestepProfile]))
        
        BackgroundLoad = [a + b for a, b in zip(BackgroundLoad, devicePower)]
    
    return BackgroundLoad

## The function SetVariableLoad creates a load profile, based on a background
## load, the devices indicated, during the time elapsed between the start and 
## end times, with the timestep indicated.

## For the function SetVariableLoad, the description of the parameters is the 
## next:
##  VariableDevices: (list) contains an array with the devices considered
##      for the background load. Each element should be a HouseholdDevice.
##  BackgroundLoad: (array) contains an array with the background load.
##  EndTime: holds the end time for the simulation, in hours.
##  TimeStep: holds the sampling time for the resulting background load profile
##      in hours. It is considered as default a 15 min (0.25 h) step.
##  StartTime: holds the end time for the simulation, in hours. It is 
##      considered that the simulation starts at 0 time.
    
## Considerations:
##  - The EndTime, TimeStep and StartTime should be consistent with the 
##    background load provided.
##  - The resulting array will have as minimum timestep the timestep used for
##    the device load profile.    
    
def SetVariableLoad(VariableDevices, BackgroundLoad, EndTime, TimeStep = 0.25, StartTime = 0):
    from numpy import append
    TotalLoad = BackgroundLoad
    DeactivatedDevice = array([0 for i in VariableDevices[0].TimestepProfile])
    
    for device in VariableDevices:
        devicePower = array([])          
    
        for instant in range(int(len(BackgroundLoad)/TimeStep/60)):
#            instant_Probability = SetProbability()
#            print('The Load probability is: ', device.HourlyProbability[instant], ', whereas the instant probability is: ', instant_Probability)
            
            if device.HourlyProbability[instant] >= SetProbability():
                devicePower = append(devicePower, array(device.TimestepProfile))                               
            else:
                devicePower = append(devicePower, DeactivatedDevice)

          
        TotalLoad = [a + b for a, b in zip(TotalLoad, devicePower)]

    
    return TotalLoad

###############################################################################
    
def Create_profile(normalize = True, yearly_consumption = 1, n_houses = 1, create_CSV = False, create_figure = False):
    
    from numpy import array
    import matplotlib.pyplot as plt

    
    ###############################################################################
    
    
    ## Some examlpe devices are loaded from the library.
    Fridge = HouseholdDevice.Fridge()
    Router = HouseholdDevice.Router()
    Dishwasher = HouseholdDevice.Dishwasher()
    Lights = HouseholdDevice.Lights()
    TV = HouseholdDevice.TV()
    Cell_charger = HouseholdDevice.Cell_charger()
    Laptop = HouseholdDevice.Laptop()
    Screen = HouseholdDevice.Screen()
    Kitchen = HouseholdDevice.Kitchen()
    Microwave = HouseholdDevice.Microwave()
    Kettle = HouseholdDevice.Kettle()
    
    ## A list of devices is created for the background profile.
    BackgroundDevices = [Fridge, Router, Dishwasher]
    
    ## A list of devices is created for the variable load profile.
    VariableDevices = [Lights, TV, Cell_charger, Laptop, Screen, Kitchen, Microwave, Kettle]
    
    ## The simulation time is set
    ndays = 365
    StartTime = 0
    EndTime = 24*ndays
    
    ## Number of apartments in the building.
    BuildingLoad_total = array([0]*60*int((EndTime-StartTime)))
    
    
    
    for i in range(n_houses):    
    
    ## For each house, the background load and the total load is calculated, and
    ## then added to the building's load.
        BackgroundLoad = []
        TotalLoad = []
        for d in range(ndays):
        
    ## Base load for the building.
            BackgroundLoad += SetBackgroundLoad(BackgroundDevices, 24)
            TotalLoad += SetVariableLoad(VariableDevices, BackgroundLoad[-24*60:], 24)
        
        BuildingLoad = [a + b for a, b in zip(BackgroundLoad, TotalLoad)]
        
    
        BuildingLoad_total = [a + b for a, b in zip(BuildingLoad_total, BuildingLoad)]
        
    ## The load is sampled every 15 min, as the load is samplead every 1 min.
    Load_15min = [BuildingLoad_total[15*i] for i in range(96*ndays)]


    if normalize:
        Base_Load = sum(Load_15min)*0.25
        normalize_ratio = yearly_consumption/Base_Load
        
        Load_15min = [i*normalize_ratio for i in Load_15min]
    
    if create_CSV:
        from csv import writer    
        
        ## A .csv file is created to store the resulting load profile.
        file = open('Load_Profile_15min_2021_2023.csv', 'w')
        writer = writer(file)
        writer.writerow(Load_15min)
        file.close()
    
    if create_figure:
    
        ### A figure shows the load profile sampled every 15 min.
        plt.figure()
        plt.plot(Load_15min, 'b')
        plt.xlim([1, len(Load_15min)])
    
    
        plt.xlabel('Timestep')
        plt.ylabel('Power [kW]')
        plt.title('Load profile of the building sampled every 15 min')
        plt.grid()
    
#        ## The total energy of the building is calculated. Note that, as the output is
#        ## on kWh, and the load profile is sampled every min, the conversion to get
#        ## kW-min to kWh is 1/60.
#        TotalEnergy = sum(BuildingLoad_total)
#        TotalEnergy /= 60                       # In kWh
#        
#        print('The total energy consumed per day by the building is: ', TotalEnergy, ' kWh')
#        print('The average energy consumed per apartment in the building is: ', TotalEnergy/n_houses, ' kWh')    
    return Load_15min