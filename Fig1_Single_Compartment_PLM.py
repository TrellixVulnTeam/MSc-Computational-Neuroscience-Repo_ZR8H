# -*- coding: utf-8 -*-
"""
Created on Fri May  8 10:28:24 2020

@author: eshor
"""

import matplotlib.pyplot as plt
import seaborn as sns
from Graphing import clcolor, kcolor, xcolor,nacolor,wcolor
import Class_Single_Compartment_PLM as PLM

#### Figure 1B

def F1b(init_cl=[1e-3,15e-3,50e-3,90e-3]):
    
    Figure_1 = plt.figure()
    
    Subplot_1 = Figure_1.add_subplot(2,1,1)
    Subplot_2 = Figure_1.add_subplot(2,1,2)
    Subplot_1.set_ylabel = '[Cl-] (mM)'
    Subplot_2.set_ylabel = 'Volume (pL)'
    Subplot_1.set_axis_off
    sns.despine
    
    for i in range(len(init_cl)):
        Fig1B_Sim = PLM.SingleCompartment_PumpLeak()
        Fig1B_Sim.Set_Concentrations(Cl=[init_cl[i],119e-3],K=[0,3.5e-3])
        Fig1B_Sim.Set_LeakConductances()
        Fig1B_Sim.Set_NaKATPase_Properties(time_off=1000000,time_on=0)
        Fig1B_Sim.Set_KCC2_Properties()
        Fig1B_Sim.Set_Timing(duration = 5000 ,time_step = 1e-3)
        Fig1B_Sim.Simulate_PLM()
        if i ==0:
            Subplot_1.plot(Fig1B_Sim.t_arr,Fig1B_Sim.cl_arr,color=clcolor,ls='solid')
            Subplot_2.plot(Fig1B_Sim.t_arr,Fig1B_Sim.volume_arr,color='k',ls='solid')
        else: 
            Subplot_1.plot(Fig1B_Sim.t_arr,Fig1B_Sim.cl_arr,color=clcolor,ls='--')
            Subplot_2.plot(Fig1B_Sim.t_arr,Fig1B_Sim.volume_arr,color='k',ls='--')
       
def F1c():
    
    
      
    Fig1C_Sim = PLM.SingleCompartment_PumpLeak()
    Fig1C_Sim.Set_Concentrations()
    Fig1C_Sim.Set_LeakConductances()
    Fig1C_Sim.Set_NaKATPase_Properties(time_off=3000,time_on=9000,switch_at_start='on')
    Fig1C_Sim.Set_KCC2_Properties()
    Fig1C_Sim.Set_Timing(duration = 12000 ,time_step = 1e-3)
    Fig1C_Sim.Simulate_PLM()
    
    Figure_1C, (ax1,ax2,ax3) = plt.subplots(nrows=3, ncols=1, sharex='all')
    
    ax1.plot(Fig1C_Sim.t_arr,Fig1C_Sim.cl_arr,color=clcolor,ls='solid',label='[Cl]')
    ax1.plot(Fig1C_Sim.t_arr,Fig1C_Sim.k_arr,color=kcolor,ls='solid', label='[K]')
    ax1.plot(Fig1C_Sim.t_arr,Fig1C_Sim.na_arr,color=nacolor,ls='solid', label = '[Na]')
    ax1.plot(Fig1C_Sim.t_arr,Fig1C_Sim.x_arr,color=xcolor,ls='solid', label ='[X]')
    ax1.set_ylabel('Concentrations(mM)')
    ax1.legend()
    

    ax2.plot(Fig1C_Sim.t_arr,Fig1C_Sim.e_cl_arr,color=clcolor,ls='solid', label='ECl')
    ax2.plot(Fig1C_Sim.t_arr,Fig1C_Sim.e_k_arr,color=kcolor,ls='solid',label='EK')
    ax2.plot(Fig1C_Sim.t_arr,Fig1C_Sim.vm_arr,color='k',ls='solid', label='Vm')
    ax2.set_ylabel('Voltages (mV)')
    ax2.legend()

    ax3.plot(Fig1C_Sim.t_arr,Fig1C_Sim.volume_arr,color='k',ls='solid', label='Volume')
    ax3.set_ylabel('Volume (pL)')
    ax3.set_xlabel('Time (ms)')
    ax3.legend()
    sns.despine()   
    
    Fig1
    

F1c()


#F1b()
