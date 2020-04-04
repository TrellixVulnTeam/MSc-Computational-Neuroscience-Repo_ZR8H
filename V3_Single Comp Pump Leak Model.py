# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:09:13 2020

@author: eshor
"""
#########################################################
####### 1) IMPORT MODULES ###########
#********************************************************
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
print('Modules Imported')
#########################################################

#########################################################
####### 2) SET CONSTANTS ##########
#********************************************************
R = 8.31466  # J/(Kelvin.mol) #Universal Gas Constant
F = 96485.0 #C/mol # Faraday Constant in Volts
T = 310.15   #Kelvin # Absolute temperature (37C)
RTF = R*T/F  # J/C
pw = 0.0015 #dm/s #osmotic permeability of biological membranes
vw = 0.018 #dm^3/mol #partial molar volume of water
print('Constants set')
print("")
#########################################################

def Set_Single_Compartment_Parameters(radius = 5e-5, length=25e-5, Capacitance=2e-4):
    w_initial = np.pi*radius**2*length #Initial volume in litres
    w = w_initial
    SA_initial = 2*np.pi*radius*length #Initial surface area in dm^2
    SA = SA_initial 
    Ar = w/SA               #Area scaling constant dm^-1
    C_m = Capacitance
    FinvCAr = F/(C_m)*Ar
    print('Initial compartment volume:', w_initial, 'Litres')
    print("")
    
def Set_Timing(t_time,):
def Set_Concentrations():
def run_Single_Compartment_plm()
    


#########################################################
######### 3) SET PARAMETERS #########
#*******************************************************
######### 3.1) DIMENSIONS AND MEMBRANE PROPERTIES
#*******************************************************


C_m = 2e-4              #(F/dm^2) #Unit of membrane capacitance

#*********************************************************
######### 3.2) IONIC CONCENTRATIONS
#********************************************************* 
#Extracellular concentrations (M)
ConcO_K = 3.5e-3             
ConcO_Cl= 119e-3      
ConcO_Na = 145e-3 
ConcO_X = -1*(ConcO_Cl - ConcO_Na - ConcO_K)*0.2
z_X = -0.85 
OsmolO = ConcO_K + ConcO_Na + ConcO_Cl + ConcO_X
#Intracellular concentrations (M)
ConcI_K = 122.873e-3
#ConcI_Cl = 5.163e-3
ConcI_Cl = 60e-3
ConcI_Na = 14.002e-3
ConcI_X = 154.962e-3
OsmolI = ConcI_K + ConcI_Na + ConcI_Cl + ConcI_X
ratio_X = 0.98 # set ratio of intracellular:extracellular impermeants
print("Extracellular concentrations: (Molar)")
print("---------------------")
print('Na: ',ConcO_Na)
print('K: ',ConcO_K)
print('Cl: ', ConcO_Cl)
print('X: ', ConcO_X)
print('Extracellular Osmolarity:', OsmolO)
print("")
print("Intracellular concentrations: (Molar)")
print("---------------------")
print('Na: ',ConcI_Na)
print('K: ',ConcI_K)
print('Cl: ', ConcI_Cl)
print('X: ', ConcI_X)
print('Intracellular Osmolarity:', OsmolI)
print("")
#****************************************************
###### 3.3) CHANNEL AND PUMP PROPERTIES ######
#****************************************************
g_Na = 2e-3/F   # S/dm2 #Na Leak conductance  (All from Kira's code)
g_K = 7e-3/F    # S/dm2 #K Leak conductance
g_Cl =2e-3/F    # S/dm2 #Cl Leak conductance
g_KCC2 = 2e-3/F #KCC2 conductance from Kira's code
g_X = 0/F       # No conductance for impermeant anions

curr = -0*5e-8; #Baseline no current injected

p_default = -1
p_effective = (10**p_default)/F

print('Pump rate:', p_effective)
#****************************************************
###### 3.4) TIMING VARIABLES
#****************************************************
t = 0               # real time in seconds
totalt = 1000        # total time
dt = 1e-3           # time step
t_on = 0            # time when the ATPase is switched on
t_off = 5000        # time when the ATPase is switched off
totalsteps = round(totalt/dt)  # total number of time steps

sw = 1              # 0 = switch off; 1 = switch on
ctr = 1     # counter for plotting 
n = 1800     # number of plot points
ts = totalt/n   # time interval for plotting

print('Simulation will run for ',totalt,' seconds') 
print('The timestep will be', dt, 'seconds')
print('There will be ', totalsteps, ' time steps')
print('The ATPase pump is switched on at', t_on, 'seconds')
print('The ATPase pump is switched off at', t_off, 'seconds')

#****************************************************
###### 3.4) ARRAY INITIALIZATION
#****************************************************
Na_Arr = []
K_Arr = []
Cl_Arr = []
X_Arr = []
Vm_Arr = []
t_Arr = []           
w_Arr = []
ECl_Arr =[]
EK_Arr = []

print('All arrays are empty to begin with')
########################################################

########################################################
# 4) SIMULATION #######
########################################################

Net_Intracellular_Charge = (ConcI_Na+ConcI_K)-ConcI_Cl+z_X*ConcI_X
Vm = FinvCAr*(Net_Intracellular_Charge)*1e3

print('Net Intracellular Charge:', Net_Intracellular_Charge)
print('F*w/Cm*SA :', FinvCAr)
print('Starting voltage of: ', round(Vm,3), 'mV')
p_default = -1
p_effective = (10**p_default)/F

print('Pump rate:', p_effective)

for i in range(1,totalsteps):  #note tit = total number of timesteps 
    
    
    #Determining switch position of ATPase
    if (t<t_off) & (t>t_on): 
        sw=1 
    else: 
        sw=0
            
    # Voltage Calculation 
    Net_Intracellular_Charge = (ConcI_Na+ConcI_K)-ConcI_Cl+z_X*ConcI_X
    Vm = FinvCAr*(Net_Intracellular_Charge)
    
    # KCC2 pump rate calculation 
    EK = RTF*np.log(ConcO_K/ConcI_K)
    ECl = -RTF*np.log(ConcO_Cl/ConcI_Cl)
    JKCC2 = g_KCC2*(EK-ECl)

    # ATPase pump rate calclation
    Jp = p_effective*(ConcI_Na/ConcO_Na)**3

    #Dummy X(impermeant anion) concentrations:
    #ConcO_X_temp = ConcI_X*(1-X_ratio)    

    #Incrementing Ion concentration    
    d_Na = -dt/Ar*(g_Na*(Vm-RTF*np.log(ConcO_Na/ConcI_Na)) + sw*3*Jp) 
    d_K = -dt/Ar*(g_K*(Vm-RTF*np.log(ConcO_K/ConcI_K))-JKCC2-sw*2*Jp) 
    d_Cl = dt/Ar*(g_Cl*(Vm+RTF*np.log(ConcO_Cl/ConcI_Cl))+JKCC2)
    d_X = -dt/Ar*z_X*(g_X*np.log(ConcO_X/ConcI_X))
    ConcI_Na += d_Na   
    ConcI_K += d_K
    ConcI_Cl += d_Cl    
     
    #Osmolarity and volume adjustments
    OsmolI = ConcI_Na+ConcI_K+ConcI_Cl+ConcI_X
    w2 = w*(OsmolI/OsmolO)
    
    #Adjusting Concentrations based on new volumes
    ConcI_Na *= w/w2
    ConcI_K *= w/w2
    ConcI_Cl *= w/w2
    ConcI_X *= w/w2
    
    #Updating Arrays and counters
    w = w2
    
    if t >= ctr*ts :
        Vm_Arr.append(Vm*1e3)
        ECl_Arr.append(ECl*1e3)
        EK_Arr.append(EK*1e3)
        K_Arr.append(ConcI_K*1e3)
        Na_Arr.append(ConcI_Na*1e3)
        Cl_Arr.append(ConcI_Cl*1e3)
        X_Arr.append(ConcI_X*1e3)
        w_Arr.append(100*(1e5)*((3/(4*np.pi))*w)**(1/3))
        t_Arr.append(t)
        ctr += 1
    
        
    t = t+dt      
print('Simulation complete')
print('Final chloride concentration:', Cl_Arr[-1], ' mM')
print('Final cell volume:', w_Arr[-1], 'L')


########################################################

########################################################
# 5) PLOTS
########################################################

"""
plt.plot(t_Arr,Na_Arr,label ="Na")
plt.plot(t_Arr,Cl_Arr,label="Cl")
plt.plot(t_Arr,K_Arr,label="K")
plt.plot(t_Arr,X_Arr,label="X")

SwitchOffX_Arr =[]
SwitchOnX_Arr= []
SwitchOnY_Arr=[]
SwitchOffY_Arr=[]
for a in range(0,150):
  SwitchOffX_Arr.append(t_off)
  SwitchOnX_Arr.append(t_on)
  SwitchOnY_Arr.append(a)
  SwitchOffY_Arr.append(a)

plt.plot(SwitchOnX_Arr,SwitchOnY_Arr,'k:')
plt.plot(SwitchOffX_Arr,SwitchOffY_Arr,'k:')
plt.title("Ionic concentration changes")
plt.xlabel("Time (s)")
plt.ylabel("Intracellular ionic concentrations (M)")
plt.legend(loc = 'upper right')
sns.despine()
plt.annotate('ATPase On', xy=(t_on, 150))
plt.annotate('ATPase Off', xy=(t_off, 150))


plt.plot(t_Arr,EK_Arr,label ="EK")
plt.plot(t_Arr,ECl_Arr,label="ECl")
plt.plot(t_Arr,Vm_Arr,label="Vm")

SwitchOffX_Arr =[]
SwitchOnX_Arr= []
SwitchOnY_Arr=[]
SwitchOffY_Arr=[]
for a in range(-120,150):
  SwitchOffX_Arr.append(t_off)
  SwitchOnX_Arr.append(t_on)
  SwitchOnY_Arr.append(a)
  SwitchOffY_Arr.append(a)

plt.plot(SwitchOnX_Arr,SwitchOnY_Arr,'k:')
plt.plot(SwitchOffX_Arr,SwitchOffY_Arr,'k:')
plt.title("Changes in cell potentials")
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.legend(loc = 'upper right')
sns.despine()
plt.annotate('ATPase On', xy=(t_on, 150))
plt.annotate('ATPase Off', xy=(t_off, 150))


plt.plot(t_Arr,w_Arr,label="Volume (pL)")


SwitchOffX_Arr =[]
SwitchOnX_Arr= []
SwitchOnY_Arr=[]
SwitchOffY_Arr=[]
for a in range(300,530):
  SwitchOffX_Arr.append(t_off)
  SwitchOnX_Arr.append(t_on)
  SwitchOnY_Arr.append(a)
  SwitchOffY_Arr.append(a)
  
plt.plot(SwitchOnX_Arr,SwitchOnY_Arr,'k:')
plt.plot(SwitchOffX_Arr,SwitchOffY_Arr,'k:')
sns.despine()

plt.annotate('ATPase On', xy=(t_on, max(w_Arr)))
plt.annotate('ATPase Off', xy=(t_off, max(w_Arr)))"""
