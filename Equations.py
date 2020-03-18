# -*- coding: utf-8 -*-
"""
This python code will be contain methods for all the electrophyiological equations in the Single Compartment Neuron))

Created on Fri Feb 21 17:03:16 2020

@author: eshor

EQUATIONS:
    
     Calc_I_Ion: Calculates the current of a particular ion
     Calc_E_Ion: Calculates the Potential of an ion at a point in time
     Calc_Vm: Calculates the membrane potential based on the charge difference approach
     Calc_Jp: Calculates the pump-rate of the Na-K-ATPase
     Calc_Jkcc2: Calculates the pump-rate of the KCC2 cotransporter
     Calc_Osmol: Calculates the internal and external osmolalities at a given time.
     Calc_dw: Calculates the change in unconstrained volume of the neuron per delta t
     Calc_Am: Calculates Am (ratio of SA to volume of compartment)
     Calc_dNa: Calculates the change of intracellular Na+ per dt
     Calc_dK: Calculates the change of intracellular K+ per dt
     Calc_dCl: Calculates the change of intracellular Cl- per dt
     Calc_Hp: Calculates Hydrostatic pressure
     Calc_dw_Constrained: Calculates the dw constrained by hydrostatic pressures
     Calc_SA: Calculates Surface Area of the cell
     Calc_w: Calculates the Volume of the cell
    """


import numpy as np

############################################################################
#******************************** CONSTANTS *******************************#
############################################################################

R = 8.31466 # J/(K.mol) #Universal Gas Constant
F = 96485.33 #C/mol # Faraday Constant in Volts
T = 310.15 #K # Absolute temperature (37C)
RTF = R*T*(1/F)

############################################################################
#******************************** METHODS *********************************#
############################################################################



def Calc_I_Ion(g_ion,V_m,E_ion):
   "Calculates the current of a particular ion"
   I_ion = g_ion*(V_m - E_ion)
   return I_ion


def Calc_E_Ion(z_ion,ConcOut_ion, ConcIn_ion):
    E_ion= RTF/z_ion*np.log(ConcOut_ion/ConcIn_ion)
    return E_ion


def Calc_Vm(ConcI_Na,ConcI_K,ConcI_Cl,z,ConcI_X,Cm,Am):
    numerator = F*(ConcI_Na+ConcI_K-ConcI_Cl+z*ConcI_X) 
    denominator = Cm*Am
    Vm = numerator/denominator *1e3 #analytical units of the formula get you to uV, this is converted to mV
    return Vm


def Calc_Jp(P,ConcI_Na,ConcO_Na):
    Jp = P * ((ConcI_Na/ConcO_Na)**3)
    return Jp


def Calc_Jkcc2(g_kcc2, E_K,E_Cl):
    Jkcc2 = g_kcc2*(E_K - E_Cl)
    return Jkcc2


def Calc_Osmol(ConcI_Na,ConcO_Na,ConcI_K,ConcO_K,ConcI_Cl,ConcO_Cl,ConcI_X,ConcO_X):
    Osmol_I = ConcI_Na + ConcI_K -ConcI_Cl-ConcI_X
    Osmol_O = ConcO_Na + ConcO_K - ConcO_Cl - ConcO_X
    return [Osmol_I,Osmol_O]


def Calc_dw(vw,pw,SA,Osmol_I, Osmol_O,dt):
    dw = vw*pw*SA*(Osmol_I-Osmol_O)*dt
    return dw


def Calc_Am(SA,w):
    Am = SA/w
    return Am


def Calc_dNa(Am,gNa,Vm,ENa,Jp,w,ConcI_Na,dt,dw):
    dNa = (-Am/F)*dt*(gNa*(Vm-ENa)+3*Jp) - (1/w)*dw*ConcI_Na
    return dNa


def Calc_dK(Am,gK,Vm,EK,Jp,Jkcc2,w,ConcI_K,dt,dw):
    dK = (-Am/F) *dt*(gK*(Vm-EK)-2*Jp-Jkcc2) - 1/w*dw*ConcI_K
    return dK


def Calc_dCl(Am,gCl,Vm,ECl,Jkcc2,w,ConcI_Cl,dt,dw):
    dCl = (Am/F) *dt*(gCl*(Vm-ECl)+Jkcc2) - 1/w*dw*ConcI_Cl
    return dCl
    

def Calc_Hp(rad_Now, rad_Rest,km):
    if rad_Now > rad_Rest:
        Hp= 4*np.pi*km*(1-rad_Rest/rad_Now)
    else:
        Hp=0
    return Hp


def Calc_dw_Constrained(vw,pw,SA,Osmol_I, Osmol_O,dt,Hp):
    dw = vw*pw*SA*(Osmol_I-Osmol_O-(Hp/(R*T)))*dt
    return dw
    

def Calc_SA(diam,length):
    rad = diam/2
    SA = 2*np.pi*rad*length #Ends have been removed
    return SA

def Calc_w(diam,length):
    rad = diam/2
    w = (np.pi*rad**2)*length #Ends have been removed
    return w

