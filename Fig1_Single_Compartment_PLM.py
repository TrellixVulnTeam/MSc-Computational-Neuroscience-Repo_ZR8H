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
"""def F1b(init_cl=[1e-3,15e-3,50e-3,90e-3]):
    
    cl_arr =[]
    volume_arr =[]
    t_arr = []
    Figure1 = plt.figure()
    for i in range(len(init_cl)):
        
        Fig1B_Sim = PLM.SingleCompartment_PumpLeak()
        Fig1B_Sim.Set_Concentrations(Cl=[init_cl[i],119e-3])
        Fig1B_Sim.Set_LeakConductances()
        Fig1B_Sim.Initialize_Arrays()
        Fig1B_Sim.Set_NaKATPase_Properties()
        Fig1B_Sim.Set_KCC2_Properties()
        Fig1B_Sim.Set_Timing()
        Fig1B_Sim.Simulate_PLM()
        #Figure1.subplot(2,1,1)
        Figure1.plot(Fig1B_Sim.t_arr,Fig1B_Sim.cl_arr,color=clcolor)
        #Figure1.subplot(2,1,2)
        #Figure1.plot(Fig1B_Sim.t_arr,Fig1B_Sim.volume_arr)
        

"""
Fig1B_Sim = PLM.SingleCompartment_PumpLeak()
Fig1B_Sim.Set_Concentrations(Cl=[15e-3,119e-3],K=[0,3.5e-3])
Fig1B_Sim.Set_LeakConductances()
Fig1B_Sim.Initialize_Arrays()
Fig1B_Sim.Set_NaKATPase_Properties(time_off=1000000,time_on=0)
Fig1B_Sim.Set_KCC2_Properties()
Fig1B_Sim.Set_Timing(duration = 5000 ,time_step = 1e-3)
Fig1B_Sim.Simulate_PLM()
plt.plot(Fig1B_Sim.t_arr,Fig1B_Sim.cl_arr,color=clcolor)
#plt.plot(Fig1B_Sim.t_arr,Fig1B_Sim.volume_arr)