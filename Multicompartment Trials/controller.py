"""
Main script to run simulation
"""
import simulator2
import compartment

# 1) CREATE COMPARTMENTS:
comp1 = compartment.Compartment("Comp1")
comp1.set_ion_properties()
comp2 = compartment.Compartment("Comp2")
comp2.set_ion_properties()
comp2.xflux()
comp3 = compartment.Compartment("Comp3")
comp3.set_ion_properties()


# 2) DEFINE SIMULATOR CLASS
sim = simulator2.simulator()
sim.add_compartment("comp1")
sim.add_compartment("comp2")
sim.add_compartment("comp3")

sim.set_external_ion_properties()
sim.set_timing(total_t=1000, time_step=0.0001)
sim.set_area_scale(constant_ar=False)
sim.set_electrodiffusion_properties(ED_on=True)
sim.set_j_atp(constant_j_atp=False)
sim.set_xflux(comps = ["comp2"],start_t=0)
sim.set_xoflux()
sim.set_zflux()



# 3) RUN SIMULATION
sim.run_simulation()

# 4) ACCESS RESULTS