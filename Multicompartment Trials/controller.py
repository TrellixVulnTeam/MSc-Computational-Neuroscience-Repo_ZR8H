"""
Main script to run simulation
"""
import time

import compartment
import simulator2
import pandas as pd
import h5py


file_name ="Experiment-B3"


# 1) DEFINE SIMULATOR CLASS AND ADD COMPARTMENTS
sim = simulator2.simulator(file_name)
comp1 = compartment.Compartment("Comp1")
comp1.set_ion_properties()
sim.add_compartment(comp1)
comp2 = compartment.Compartment("Comp2")
comp2.set_ion_properties()
sim.add_compartment(comp2)
comp3 = compartment.Compartment("Comp3")
comp3.set_ion_properties()
sim.add_compartment(comp3)
comp4 = compartment.Compartment("Comp4")
comp4.set_ion_properties()
sim.add_compartment(comp4)
comp5 = compartment.Compartment("Comp5")
comp5.set_ion_properties()
sim.add_compartment(comp5)


#sim.add_default_multicompartment(number_of_comps=9)

# 2) SET SIMULATION SETTINGS

sim.set_electrodiffusion_properties(ED_on=True)

sim.set_external_ion_properties()
sim.set_j_atp(constant_j_atp=False)
sim.set_area_scale(constant_ar=False)
total_t = 6*60
time_step = 1e-6
sim.set_timing(total_t=total_t, time_step=time_step, intervals=1000)

sim.set_xflux(comps=["Comp2"], flux_type="static", start_t=60, end_t=180, x_conc=1e-3, z=-2.0, flux_rate=5*1e-3/60)
#sim.set_xflux(comps=["Comp3"], flux_type="static", start_t=100, end_t=500, x_conc=1e-3, z=-1.0, flux_rate=0.4*1e-3/60)
# sim.set_zflux()
##sim.set_xoflux()
#run_t_arr = np.arange[0:total_t:time_step]

sim.run_simulation()
print("fin")
# 4) RUN SIMULATION


# 5) ACCESS RESULTS

