# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 18:10:40 2021

Controls the processes of the multicompartmental model

@author: eshor
"""
import numpy as np
from matplotlib import pyplot as plt
import common
import electrodiffusion
import constants
from compartment import Compartment


total_t = 120 #total time in seconds
run_t = 0 # running time, starting at 0 seconds
dt = 1e-3 # time step in seconds
t_arr = []
run_t=0

Comp_a = Compartment("comp_a")
Comp_a.set_ion_properties()
Comp_b = Compartment("comp_b")
Comp_b.set_ion_properties()


while run_t < total_t:
         
    Comp_a.step(dt)
    Comp_a.update_volumes()
    Comp_a.update_arrays()
   # Comp_b.step(dt)
   # Comp_b.update_volumes()
   # Comp_b.update_arrays() }}}
    t_arr.append(run_t)
    run_t += dt
    
    
    
plt.plot(t_arr, Comp_a.na_arr)