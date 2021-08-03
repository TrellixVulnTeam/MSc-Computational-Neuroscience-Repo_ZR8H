"""
Version of the controller class to be able to start a simulation from where an HDF file started off
"""
import compartment
import simulator3
import h5py
# Base file is the original file where the simulation was run

base_file_name = "Experiment-D6"
new_file_name = "Test-inhib3"

sim = simulator3.simulator(new_file_name)

with h5py.File(base_file_name, mode='r') as hdf:
    C = hdf.get('COMPARTMENTS')
    C_group_arr = []
    comp_names_arr = list(C.keys())

    T = hdf.get('TIMING')
    total_t = T.get('TOTAL_T')[()]
    intervals = T.get('INTERVALS')[()]
    dt = T.get("DT")[()]
    total_steps = total_t / dt
    interval_step = total_steps / intervals
    interval_arr = [round(interval_step * i) for i in range(intervals)]

    for i in comp_names_arr:
        temp_arr =[]
        comp_name = i
        current_comp = C.get(comp_name)
        last_time_point = str(interval_arr[-1])
        last_dataset = current_comp.get(last_time_point)
        last_dataset = list(last_dataset)

        radius = last_dataset[1]
        length = last_dataset[2]
        na_i = last_dataset[4]
        k_i = last_dataset[5]
        cl_i = last_dataset[6]
        x_i = last_dataset[7]
        z_i = last_dataset[8]

        comp = compartment.Compartment(comp_name,radius,length)
        comp.set_ion_properties(na_i,k_i,cl_i,x_i,z_i,osmol_neutral_start=False)
        sim.add_compartment(comp)
        

sim.set_electrodiffusion_properties(ED_on=True)

sim.set_external_ion_properties()
sim.set_j_atp(constant_j_atp=True)
sim.set_area_scale(constant_ar=False)
total_t = 0.2
time_step = 1e-6
sim.set_timing(total_t=total_t, time_step=time_step, intervals=1000)
sim.add_synapse("Comp8", "Inhibitory", start_t=0.05, duration=2e-3, max_neurotransmitter= 50e-3)

sim.run_simulation()
print("fin")

