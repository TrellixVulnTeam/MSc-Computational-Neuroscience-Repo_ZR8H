# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 18:10:40 2021

Controls the processes of the multi-compartmental model

@author: E Shorer
"""

import common
import constants
import compartment
import electrodiffusion
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


dt = 1e-3  # 1ms time steps
total_t = 1200  # s
run_t = 0  # current simulation timing
t_arr = [0]


comp_1 = compartment.Compartment("comp_1")
comp_1.set_ion_properties()
comp_2 = compartment.Compartment("comp_2")
comp_2.set_ion_properties()
ed_1_2 = electrodiffusion.Electrodiffusion(comp_1, comp_2)

while run_t < total_t:
    comp_1.step(dt)
    comp_2.step(dt)
    comp_1_ed_dict = comp_1.get_ed_dict()
    comp_2_ed_dict = comp_2.get_ed_dict()
    ed_conc_changes = ed_1_2.calc_ed(dt, comp_1_ed_dict, comp_2_ed_dict)
    comp_1.ed_update(ed_conc_changes)
    for j in ed_conc_changes:
        ed_conc_changes[j] *= -1
    comp_2.ed_update(ed_conc_changes)
    comp_1.update_volumes()
    comp_1.update_arrays()
    run_t += dt
    t_arr.append(run_t)

plt.plot(t_arr,comp_1.k_arr)
