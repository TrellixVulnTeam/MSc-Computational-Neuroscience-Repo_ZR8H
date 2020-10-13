# -*- coding: utf-8 -*-
"""
Created on Tue May 19 16:14:15 2020


@author: eshor
"""


import matplotlib.pyplot as plt
import seaborn as sns
from Graphing import clcolor, kcolor, xcolor,nacolor
#import Class_Single_Compartment_PLM as PLM
import Class_Single_Compartment_Analytical as analytical_PLM

### FIGURE 2 ###

#Looping through conductances and getting the voltages and volumes, using the analytical solution


def F2a():

    vm_arr =[]
    ek_arr =[]
    ecl_arr =[]
    volume_arr=[]   
    g_arr = []

    for g in range(1,100):
        
        Fig2a_Sim = analytical_PLM.SingleCompartment_PumpLeak_Analytical()
        gk =g*(1e-4)
        
        Fig2a_Sim.Set_External_Concentrations()
        Fig2a_Sim.Set_Leak_Conductances(g_k= gk)
        Fig2a_Sim.Set_NaKATPase_Properties()
        Fig2a_Sim.Set_KCC2_Properties()
        Fig2a_Sim.Calc_Theta()
        Fig2a_Sim.Calc_Vm()
        Fig2a_Sim.Calc_EK()
        Fig2a_Sim.Calc_ECl()
        
            
        vm_arr.append(Fig2a_Sim.vm*1e3)
        ek_arr.append(Fig2a_Sim.Ek*1e3)
        ecl_arr.append(Fig2a_Sim.ECl*1e3)
        g_arr.append(g)
        
        
    Figure_2A, ax1 = plt.subplots(1,1)
    ax1.plot(g_arr,vm_arr, color="k")
    ax1.plot(g_arr,ek_arr, color = kcolor)
    ax1.plot(g_arr,ecl_arr, color = clcolor)
    print(vm_arr)
        
    
F2a()

