# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 16:18:48 2020

This is a Simplified Version of the Pump Leak Model found in Kira's code.
This version is aligned with Alan Kay's code. 

I will try make this code somewhat modular and increase the complexity in a step wise fashion. 

Step 1) Create the compartment and control for ionic fluxes via diffusion - done
Step 2) Add the ATP-ase pump - done
Step 3) Add the KCC2 pump - done



@author: eshor
"""


import numpy as np
import matplotlib.pyplot as plt
import Equations as eq

############################################################################
#********************************PARAMETERS *******************************#
############################################################################

R = 8.31466 # J/(Kelvin.mol) #Universal Gas Constant
F = 96485.33 #C/mol # Faraday Constant in Volts
T = 310.15 #Kelvin # Absolute temperature (37C)
RTF = R*T/F

#########  MEMBRANE and CELLULAR properties (capacitance, permeabilities, volumes) """

rad = 5*1e-5            #Radius of the single compartment (um converted to dm)
w = (4/3)*np.pi*rad**3  #L #Volume of a single compartment (Sphere)
SA = 4*np.pi*rad**2     #dm^2 #Surface area of a sphere 
C_m = 1e-4              #(F/dm^2) #Unit of membrane capacitance
FinvCSA = F/(C_m*SA)
    
######### IONIC properties (concentrations and osmolalities)

OsmolO = 300e-3     #M #Fixed osmolarity

ConcO_K = 3e-3      #M 
ConcO_X = 0         #M #No anions outside cell
zO_X = 0            #charge of extracellular anions
ConcO_Cl=(OsmolO-ConcO_X+zO_X*ConcO_X)/2 # calc clo from charge and osmotic balance 
ConcO_Na = ConcO_Cl-ConcO_K-ConcO_X*zO_X # calc nao 
ConcI_X = 50*1e-3   #Concentration of intracellular anion (M)
ConcI_Cl = (OsmolO - 2*ConcI_X)/2
ConcI_Na = 0.8*(OsmolO-ConcI_Cl-ConcI_X)
ConcI_K = 0.2*(OsmolO-ConcI_Cl-ConcI_X)

######### CHANNEL properties (conductances, pump rate constant) 

g_Na = 0.01*0.1/F # S/dm2 #Na Leak conductance ### Not too sure why I need to / by F (according to Kira and Alan)
g_K = 0.3*0.1/F     # S/dm2 #K Leak conductance
g_Cl = 0.2*0.1/F   # S/dm2 #Cl Leak conductance

curr = -0*5e-8; #Baseline no current injected
p = 0.5e-4/F # Alan's pump rate (C/dm^2)
g_KCC2 = 2e-3/F #KCC2 conductance from Kira's code

######### TIMING variables


#time(1) = 0         # time arr starts at 0
t = 1               # real time in seconds
totalt = 10000          # total time
dt = 1e-3           # time step
t_on = 0            # time when the ATPase is switched on
t_off = 5000        # time when the ATPase is switched off
totalsteps = round(totalt/dt)  # total number of time steps
sw = 0              # 0 = switch off; 1 = switch on
ctr = 1     # counter for plotting 
n = 200     # number of plot points
ts = totalt/n   # time interval for plotting

########  ARRAY INITIALIZATION:
    
Na_Arr = []
K_Arr = []
Cl_Arr = []
X_Arr = []
Vm_Arr = []
t_Arr = []           
w_Arr = []


############################################################################
#****************************** SIMULATIONS *******************************#
############################################################################

#Starting voltage
Vm = FinvCSA*w*(ConcI_Na+ConcI_K-ConcI_Cl-ConcI_X) 


for i in range(1,totalsteps):  #note tit = total number of timesteps
     
    #Determining switch position of ATPase
    if (t<t_off) & (t>t_on): 
        sw=1 
    else: 
        sw=0
            
                    
            
    # Voltage Calculation 
    Vm = FinvCSA*w*(ConcI_Na+ConcI_K-ConcI_Cl-ConcI_X) 
    EK = RTF*np.log(ConcO_K/ConcI_K)
    ECl = RTF*np.log(ConcO_Cl/ConcI_Cl)
    JKCC2 = g_KCC2*(EK-ECl)
        
    #Incrementing Ion concentration    
    d_Na = -dt*SA*(1/w)*(g_Na*(Vm-RTF*np.log(ConcO_Na/ConcI_Na)) + sw*3*p) 
    d_K = -dt*SA*(1/w)*(g_K*(Vm-RTF*np.log(ConcO_K/ConcI_K)-JKCC2)+sw*2*p+sw*curr) 
    d_Cl = dt*SA*(1/w)*(g_Cl*(Vm+RTF*np.log(ConcO_Cl/ConcI_Cl))+JKCC2)
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
        K_Arr.append(ConcI_K*1e3)
        Na_Arr.append(ConcI_Na*1e3)
        Cl_Arr.append(ConcI_Cl*1e3)
        X_Arr.append(ConcI_X*1e3)
        w_Arr.append(100*(1e5)*((3/(4*np.pi))*w)**(1/3))
        t_Arr.append(t)
        ctr += 1
    
        
    t = t+dt      
    
####### End of simulation

####### PLOTS #######

plt.plot(t_Arr,Na_Arr,t_Arr,Cl_Arr,t_Arr,K_Arr,t_Arr,X_Arr)



    