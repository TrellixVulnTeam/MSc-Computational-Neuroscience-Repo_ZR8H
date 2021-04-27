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
comp3 = compartment.Compartment("Comp3")
comp3.set_ion_properties()

# 2) DEFINE SIMULATOR CLASS AND ADD COMPARTMENTS
sim = simulator2.simulator()
sim.add_compartment(comp1)
sim.add_compartment(comp2)
sim.add_compartment(comp3)

# 3) SET SIMULATION SETTINGS
sim.set_electrodiffusion_properties(ED_on=True)
sim.set_external_ion_properties()
sim.set_j_atp(constant_j_atp=False)
sim.set_area_scale(constant_ar=False)
sim.set_timing(total_t=100, time_step=0.001)

sim.set_xflux(comps=["Comp2"], type="dynamic", start_t=50, end_t=80, x_conc=2e-3, z=-0.85)
# sim.set_zflux()
##sim.set_xoflux()

# 4) RUN SIMULATION
sim.run_simulation()



# 5) ACCESS RESULTS

