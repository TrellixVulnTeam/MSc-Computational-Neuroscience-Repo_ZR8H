"""" SINGLE COMPARTMENT PUMP-LEAK NEURON
 Author: Eran Shorer
 Last Update on: 20/02/2020 @ 19.47

Questions for Kira:
1) Why is R (Universal gas constant) in her single compartment = 26.725*1e-3?

2) What is going on with the osmolarity calculations
    
3) Should I rump the pumps before or after the leak channel?
    
    
    ** Currently stuck at why the Vm isn't updating properly

 """


import numpy as np
import matplotlib.pyplot as plt
#from Equations import Calc_Am, Calc_Vm, Calc_DeltaConc, Calc_I_Ion, Calc_E_Ion, Calc_Jp, Calc_Jkcc2, Calc_Ol
import Equations as eq

############################################################################
#******************************** CONSTANTS *******************************#
############################################################################

R = 8.31466 # J/(K.mol) #Universal Gas Constant
F = 96485.33 #C/mol # Faraday Constant in Volts
T = 310.15 #K # Absolute temperature (37C)
RTF = R*T*(1/F)

######### CHANNEL properties (conductances, pump rate constant) 

g_Na = 2e-3 # Siemens #Na Leak conductance
g_K = 7e-3 # S #K Leak conductance
g_Cl = 2e-3 # S #Cl Leak conductance
g_KCC2 = 2e-3 # S #KCC2 conductance

p = 0.1 # C/(dm^2.s) #default pump rate constant

######### IONIC properties (concentrations and osmolalities)

ConcI_K = 122.9*1e-3 #M #Intracellular K+ concentration
ConcI_Na = 14*1e-3 #M #Intracellular Na+ concentration
ConcI_Cl = 5.2*1e-3 #M #Intracellular Cl- concentration
ConcI_X = 154.9*1e-3 #M #Intracellular impermeant anion (X) concentration

ConcO_Na = 145*1e-3 #M #Extracellular Na+ concentration
ConcO_K = 3.5*1e-3 #M #Extracellular K+ concentration
ConcO_Cl = 119*1e-3 #M #Extracellular Cl- concentration
ConcO_X = 29.5*1e-3 #M #Extracellular impermeant anion (X) concentration
testOsm = ConcO_Na+ConcO_K-ConcO_Cl-ConcO_Cl

z = -0.85 #Mean charge of intracellular and extracellular anions X

#########  MEMBRANE and CELLULAR properties (capacitance, permeabilities, volumes) """
C_m = 2*10e-4 # F/dm^2 #Unit of membrane capacitance
k_m = 25 #N/dm #Variable for membrane tension
w = 2*1e-12 #L #Volume of a single compartment
diam = 10*1e-5 #Diameter of the single compartment (um converted to dm)
length = 25*1e-5 #Length of the single compartment (um converted to dm)
v_H2O = 0.018 # dm^3 / mol # Partial molar volume of water
p_H2O = 0.0015 # dm/s # Osmotic permeability

v_max =5*1e-3 #M/s Vmax (Raimondo 2012)


dt = 1 #s #Time step
dv = 0 #mV #Voltage changes


E_K = 0
E_Na =0
E_Cl =0


t_Arr =[]

############################################################################
#****************************** SIMULATIONS *******************************#
############################################################################


### REPLICATING CHLORIDE LEAK - KIRA FIGURE 1B
## SIMULATION RUNNING FOR 20 MINS WITH 0.5s TIME STEP
ConcI_Cl = 50*1e-3 #M #Intracellular Cl- concentration
ConcO_Cl = 119*1e-3 #M #Extracellular Cl- concentration
dt = 0.05 #s #Time step
V_m = -72.6*1e-3 #V
V_m_Arr =[]
ConcI_Cl_arr =[]
ConcI_Na_arr =[]
ConcI_K_arr =[]

E_Cl_arr =[]

Delta_Cl_arr =[]
Delta_Na_arr =[]
Delta_K_arr =[]
w_Arr =[]


Jp_arr =[]
Jkcc2_arr = []

for t in np.arange(1,1200,dt):
    
    #Calculating the SA, Am (ratio of SA:w), Vm,dw(delta volume):
    SA = eq.Calc_SA(diam, length)
    A_m = eq.Calc_Am(SA,w)
    V_m = eq.Calc_Vm(ConcI_Na,ConcI_K,ConcI_Cl,z,ConcI_X,C_m,A_m)
    V_m_Arr.append(V_m)
    
    #Calculating internal and external osmolalities and change in cell volume
    #OsmolI,OsmolO = eq.Calc_Osmol(ConcI_Na,ConcO_Na,ConcI_K,ConcO_K,ConcI_Cl,ConcO_Cl,ConcI_X,ConcO_X)
    OsmolI = ConcI_Na+ConcO_K-ConcI_Cl-ConcI_X
    OsmolO = ConcO_Na+ConcO_K-ConcO_Cl-ConcO_X
    dw = eq.Calc_dw(v_H2O,p_H2O,SA,OsmolI,OsmolO,dt)
    w = w+dw
    w_Arr.append(w*1e12)
    
    #Calculating pump rates of ATPase and KCC2
    Jp = eq.Calc_Jp(p,ConcI_Na,ConcO_Na) ## XXX pump rate seems odd
    Jkcc2 = eq.Calc_Jkcc2(g_KCC2,E_K,E_Cl)
    
    #Calculating Cl changes
    E_Cl = eq.Calc_E_Ion(-1,ConcO_Cl,ConcI_Cl) ## Seems wrong
    E_Cl_arr.append(E_Cl*1e3)
    I_Cl = eq.Calc_I_Ion(g_Cl, V_m, E_Cl)
    dCl = eq.Calc_dCl(A_m, g_Cl, V_m, E_Cl, Jkcc2, w, ConcI_Cl, dt, dw)
    Delta_Cl_arr.append(dCl*1e3)
    ConcI_Cl = ConcI_Cl + dCl*1e-3
    ConcI_Cl_arr.append(ConcI_Cl*1e3)
    
    #Calculating Na changes
    E_Na = eq.Calc_E_Ion(1,ConcO_Na,ConcI_Na)
    I_Na = eq.Calc_I_Ion(g_Na, V_m, E_Na)
    dNa = eq.Calc_dNa(A_m, g_Na, V_m, E_Na, Jp, w, ConcI_Na, dt, dw)
    Delta_Na_arr.append(dNa*1e3)
    ConcI_Na = ConcI_Na + dNa*1e-3
    ConcI_Na_arr.append(ConcI_Na*1e3)
    
    #Calculating K+ changes
    E_K = eq.Calc_E_Ion(1,ConcO_K,ConcI_K)
    I_K = eq.Calc_I_Ion(g_K, V_m, E_K)
    dK = eq.Calc_dK(A_m, g_K, V_m, E_K, Jp, Jkcc2, w, ConcI_K, dt, dw)
    Delta_K_arr.append(dK*1e3)
    ConcI_K = ConcI_K + dK*1e-3
    ConcI_K_arr.append(ConcI_K*1e3)
    
    
    t_Arr.append(t)
    
    
plt.plot(t_Arr,ConcI_Cl_arr)
    
    
    
    
    


