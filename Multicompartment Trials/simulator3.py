"""

Control the functioning of the Neural Physiological Emulator


"""

import time

import h5py
import numpy as np
import pandas as pd

import compartment
import electrodiffusion

from common import F, gna,gcl,gk,gx,g_na_k_atpase,g_kcc2,cm


class simulator:

    def __init__(self, file_name=""):
        """ Compartments array needs to be in the format of compartment class"""

        self.file_name = file_name
        self.file_name = "\ " + file_name

        try:
            with h5py.File(self.file_name, mode='w') as self.hdf:
                self.hdf.create_group('COMPARTMENTS')
                self.hdf.create_group('ELECTRODIFFUSION')
                self.hdf.create_group("TIMING")
                self.hdf.create_group("X-FLUX-SETTINGS")
                self.hdf.create_group("SYNAPSE-SETTINGS")
                print("simulation file ('" + file_name + "') created in base directory")

        except:
            raise Exception("File not created")

        self.num_comps = 0
        self.one_percent_t = 0
        self.interval_num = 1
        self.steps = 0
        self.ED_on = True
        self.constant_j_atp, self.constant_ar = False, False
        self.comp_arr, self.ed_arr = [], []
        self.external_xflux_setup, self.xflux_setup, self.zflux_setup = True, True, True
        self.na_o, self.k_o, self.cl_o, self.x_o, self.z_o, self.osm_o = 0, 0, 0, 0, 0, 0
        self.p = 0
        self.start_t, self.end_t, self.run_t, self.total_t, self.dt = 0, 0, 0, 0, 0
        self.xflux_dict = {}
        self.xflux_count = 0
        self.xoflux_switch = False
        self.xoflux_params = {"start_t": 0, "end_t": 0, "xo_conc": 0, "zo": 0}
        self.xoflux_setup = True
        self.xo_start, self.cl_o_start, self.d_xoflux, self.xo_final, self.xo_flux, self.t_xoflux = 0, 0, 0, 0, 0, 0
        self.xoflux_points, self.dt_xoflux, self.xo_alpha, self.xo_beta = 0, 0, 0, 0
        self.synapse_dict = {}


    def add_compartment(self, comp=compartment):
        """Every compartment created needs to be added to the simulator"""

        new_comp = comp.get_array(time=0)
        with h5py.File(self.file_name, mode='a') as self.hdf:
            group = self.hdf.get('COMPARTMENTS')
            subgroup = group.create_group(name=comp.name)
            subgroup.create_dataset(name='0', data=new_comp)

        self.num_comps += 1
        self.comp_arr.append(comp)

    def add_default_multicompartment(self, number_of_comps=9):
        """Sets the simulation to run with the default multicompartment model -- 9 compartments + 1 soma"""

        for i in range(number_of_comps):
            comp = compartment.Compartment("Comp" + str(i + 1), radius=0.5e-5, length=10e-5)
            comp.set_ion_properties()
            self.add_compartment(comp)

        soma = compartment.Compartment("0_Soma", radius=1e-5, length=20e-5)
        soma.set_ion_properties(na_i=0.013995241563512785,k_i=0.12286753014443351,cl_i=0.005171468255812758,
                                x_i=0.15496634531836323)
        self.add_compartment(soma)


    def get_starting_df(self):
        """Function which when called will return a dataframe of the starting values for each compartment"""

        with h5py.File(self.file_name, mode='r') as hdf:

            C = hdf.get('COMPARTMENTS')
            C_group_arr = []
            comp_names_arr = list(C.keys())
            master_arr = []

            ##### LOADING COMPARTMENT DATA
            for e in range(len(comp_names_arr)):
                C_group = C.get(comp_names_arr[e])
                C_group_arr.append(C_group)
                data_arr_2 = []
                dataset = C_group.get('0')
                data_arr = []
                for d in range(len(list(dataset))):
                    data_arr.append(dataset[d])
                data_arr_2.append(data_arr)
                master_arr.append(data_arr_2)

        df_start_data = [master_arr[i][0][1:9] for i in range(len(comp_names_arr))]
        df_start = pd.DataFrame(data=df_start_data, index=comp_names_arr)
        df_start.columns = ['Radius', 'Length', 'Volume', 'Na_i', 'K_i', 'Cl_i', 'X_i', 'z_i']
        return df_start

    def set_electrodiffusion_properties(self, ED_on=True):
        self.ED_on = ED_on
        with h5py.File(self.file_name, mode='a') as self.hdf:
            comp_group = self.hdf.get('COMPARTMENTS')
            ed_group = self.hdf.get('ELECTRODIFFUSION')

            for e in range(self.num_comps - 1):
                name = self.comp_arr[e].name + ' <- ' + self.comp_arr[e + 1].name
                ed_group.create_group(name)
                comp_a = comp_group.get(self.comp_arr[e].name)

                data_a = comp_a.get('0')
                length_a = data_a[2]
                comp_b = comp_group.get(self.comp_arr[e + 1].name)
                data_b = comp_b.get('0')
                length_b = data_b[2]

                ed = electrodiffusion.Electrodiffusion(comp_a_name=self.comp_arr[e].name, comp_a_length=length_a,
                                                       comp_b_name=self.comp_arr[e + 1].name, comp_b_length=length_b)
                self.ed_arr.append(ed)

            self.hdf.close()

    def set_external_ion_properties(self, na_o=145e-3, k_o=3.5e-3, cl_o=119e-3, x_o=29.5e-3, z_o=-0.85):
        """
        Capacity to change the extracellular bath properties before the simulation
        """
        self.na_o, self.k_o, self.cl_o, self.x_o, self.z_o = na_o, k_o, cl_o, x_o, z_o
        self.x_o = -1 * (self.cl_o - self.na_o - self.k_o)
        self.osm_o = self.x_o + self.na_o + self.cl_o + self.k_o

    def set_j_atp(self, constant_j_atp=False, p=-1):
        self.constant_j_atp = constant_j_atp
        self.p = (10 ** p) / F

    def set_area_scale(self, constant_ar=False):
        self.constant_ar = constant_ar

    def set_timing(self, total_t, time_step, intervals=1000):
        self.total_t, self.dt = total_t, time_step
        with h5py.File(self.file_name, mode='a') as self.hdf:
            timing = self.hdf.get("TIMING")
            timing.create_dataset("DT", data=time_step)
            timing.create_dataset("TOTAL_T", data=total_t)
            timing.create_dataset("INTERVALS", data=intervals)

        self.total_steps = self.total_t / self.dt

        self.output_intervals = (0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75, 1)
        self.output_arr = [round(self.output_intervals[a] * self.total_steps, 0) for a in
                           range(len(self.output_intervals))]
        self.output_arr = tuple(self.output_arr)

        self.interval_step = self.total_steps / intervals
        self.interval_arr = [round(self.interval_step * i) for i in range(intervals)]
        self.interval_arr = tuple(self.interval_arr)

    # print(self.interval_arr)

    def set_xflux(self, all_comps=False, comps=None, flux_type='static', start_t=0, end_t=0, x_conc=1e-3,
                  flux_rate=1e-5,
                  z=-0.85):

        xflux_data_arr = []  # array which will be sent to the HDF5 file
        xflux_names_arr = []
        if all_comps == False:
            xflux_data_arr.append(0)  # value 0 means that all compartments is not selected
        else:
            xflux_data_arr.append(1)

        if comps is None or comps == []:
            xflux_data_arr.append(-1)  # value -1 means no compartment has been specified
        else:
            for i in range(len(self.comp_arr)):
                if comps[0] == self.comp_arr[i].name:
                    xflux_data_arr.append(i)

        if flux_type == 'static':
            xflux_data_arr.append(0)
            xflux_data_arr.append(flux_rate)
        elif flux_type == 'dynamic':
            xflux_data_arr.append(1)
            xflux_data_arr.append(x_conc)

        xflux_data_arr.append(z)
        xflux_data_arr.append(start_t)
        xflux_data_arr.append(end_t)

        for i in range(len(comps)):
            for j in range(len(self.comp_arr)):
                if comps[i] == self.comp_arr[j].name or all_comps:
                    self.comp_arr[j].xflux_switch = True
                    self.comp_arr[j].xflux_params["type"] = flux_type  # whether xflux is static or dynamic
                    self.comp_arr[j].xflux_params["start_t"] = start_t
                    self.comp_arr[j].xflux_params["end_t"] = end_t
                    self.comp_arr[j].xflux_params["x_conc"] = x_conc
                    self.comp_arr[j].xflux_params["z"] = z
                    self.comp_arr[j].xflux_params["flux_rate"] = flux_rate

        xflux_names_arr.append("X-FLUX-" + str(self.xflux_count))  # names of the xflux
        self.xflux_count +=1
        with h5py.File(self.file_name, mode='a') as self.hdf:
            xflux_group = self.hdf.get("X-FLUX-SETTINGS")
            xflux_group.create_dataset(name=xflux_names_arr[-1], data=xflux_data_arr)

    def set_zflux(self, all_comps=False, comps=None, start_t=0, end_t=0, z_end=-1):
        """
        : param comps takes a list of compartment names
        """
        for i in range(len(comps)):
            for j in range(len(self.comp_arr)):
                if comps[i] == self.comp_arr[j].name or all_comps:
                    self.comp_arr[j].zflux_switch = True
                    self.comp_arr[j].zflux_params = {"start_t": start_t, "end_t": end_t, "z": z_end}

    def set_xoflux(self, start_t=0, end_t=50, xo_conc=1e-3, z=-0.85):
        """
        CHANGE THE FLUX OF IMPERMEANTS OUTSIDE THE COMPARTMENT
        """
        self.xoflux_switch = True
        self.xoflux_params = {"start_t": start_t, "end_t": end_t, "xo_conc": xo_conc, "zo": z}

    def gen_comps(self, comp=[compartment]):
        for i in range(len(comp)):
            yield (comp[i])

    def gen_run_t(self, run_time_arr=[]):
        for i in run_time_arr:
            yield run_time_arr[i]

    def xoflux(self):

        if self.xoflux_setup:
            # starting values for external x flux
            self.xo_start = self.x_o
            self.cl_o_start = self.cl_o
            self.xo_final = self.x_o + self.xoflux_params["xo_conc"]
            self.xoflux_points = (self.xoflux_params["end_t"] - self.xoflux_params["start_t"]) * (1 / self.dt)
            self.dt_xoflux = 4 / self.xoflux_points
            self.xo_alpha = 1
            self.xo_beta = -1
            self.xoflux_setup = False

        elif self.x_o <= self.xo_final:

            self.d_xoflux = self.xo_alpha - np.e ** (self.xo_beta * self.t_xoflux)
            self.xo_flux = self.d_xoflux * self.xoflux_params["xo_conc"]
            self.x_o = self.xo_start + self.xo_flux
            self.cl_o = self.cl_o_start - self.xo_flux  # balancing the charges added externally
            self.t_xoflux += self.dt_xoflux

        else:
            return

    def add_synapse(self, comp_name='', synapse_type='Inhibitory', start_t=0, duration=2 * 1e-3,
                    max_neurotransmitter=1e-3):
        """

        @param synapse_type: either 'Inhibitory' (GABAergic) or 'Excitatory' (Glutamatergic)
        @param comp: compartment name on which to synapse onto
        @param start_t: start time for synaptic input
        @param duration: duration of synaptic input
        @param max_neurotransmitter: max neurotransmitter concentration
        @return:
        """
        self.syn_dict ={}


        for i in range(len(self.comp_arr)):
            if comp_name == self.comp_arr[i].name:
                comp_num = i
                self.syn_dict["compartment"] = comp_num

        if synapse_type == "Inhibitory":
            self.syn_dict["synapse_type"] = 0
        elif synapse_type == "Excitatory":
            self.syn_dict["synapse_type"] = 1

        self.syn_dict["start_t"] = start_t
        self.syn_dict["duration"] = duration
        self.syn_dict["end_t"] = start_t + duration
        self.syn_dict["max_neurotransmitter_conc"] = max_neurotransmitter

        self.comp_arr[comp_num].set_synapse(synapse_type, start_t, duration, max_neurotransmitter)

        with h5py.File(self.file_name, mode='a') as self.hdf:
            synapse_name = "SYNAPSE-" + self.comp_arr[comp_num].name
            syn_data_arr = list(self.syn_dict.values())
            synapse_group = self.hdf.get("SYNAPSE-SETTINGS")
            synapse_group.create_dataset(name=synapse_name, data=syn_data_arr)

        return

    def run_simulation(self):


        self.start_t = time.time()
        for i in range(len(self.comp_arr)):
            self.comp_arr[i].dt = self.dt

        while self.run_t < self.total_t:

            if self.ED_on:

                for a in self.gen_comps(self.comp_arr):

                    a.step(dt=self.dt,
                           na_o=self.na_o, k_o=self.k_o, cl_o=self.cl_o,
                           constant_j_atp=self.constant_j_atp,
                           p=self.p, p_kcc2=2e-3 / F)

                    # step for each compartment

                    if a.xflux_switch and \
                            (a.xflux_params["start_t"] <= self.run_t <= a.xflux_params["end_t"]):
                        a.x_flux()

                    if a.zflux_switch and \
                            (a.zflux_params["start_t"] <= self.run_t <= a.zflux_params[
                                "end_t"]):
                        a.z_flux()

                    if self.xoflux_switch and \
                            self.xoflux_params["start_t"] <= self.run_t <= self.xoflux_params["end_t"]:
                        self.xoflux()

                    if a.synapse_on:
                        if self.run_t >= self.syn_dict['start_t'] and self.run_t <= self.syn_dict['end_t']:
                            a.synapse_step(run_t=self.run_t)
                    # electrodiffusion dictionary for each compartment

                for b in range(len(self.ed_arr)):
                    ed_conc_changes = self.ed_arr[b].calc_ed(self.dt, self.comp_arr[b].w,
                                                             self.comp_arr[b].get_ed_dict(),
                                                             self.comp_arr[b + 1].get_ed_dict())
                    self.comp_arr[b].ed_update(ed_conc_changes, "positive")
                    self.comp_arr[b + 1].ed_update(ed_conc_changes, "negative")

                    # appending the electrodiffusion concentrations for each compartment


                for d in self.gen_comps(self.comp_arr):
                    d.update_volumes(self.dt, self.osm_o,
                                     self.constant_ar)  # updates of the volumes, arrays, and dataframe for each compartment

                for f in range(len(self.output_arr) - 1):
                    if self.steps == self.output_arr[f]:
                        if f == 2:
                            self.one_percent_t = time.time() - self.start_t
                            self.hundred_percent_t = self.one_percent_t * 100
                            print(str(self.output_intervals[f] * 100) + " % complete in " + str(
                                round(self.one_percent_t, 2)) + " s")
                            print("Estimated time to complete :" + str(
                                round(self.hundred_percent_t / 60, 2)) + " minutes")
                        else:
                            print(str(self.output_intervals[f] * 100) + " % complete in " + str(
                                round(time.time() - self.start_t, 2)) + " s")

                if self.interval_num < len(self.interval_arr):
                    if self.steps == self.interval_arr[self.interval_num]:
                        self.interval_num += 1
                        self.save_to_file()

            self.steps += 1

            self.run_t += self.dt

        print("100.0 % complete in " + str(
            round(time.time() - self.start_t, 2)) + " s")
        self.end_t = time.time()

    """elif not self.ED_on:  # if you want to run with normal diffusion not ED
            for a in range(len(self.comp_arr)):
                self.comp_arr[a].step()
                self.comp_arr[a].x_flux()
                self.comp_arr[a].update_volumes(ar_constant)  # updates of the volumes, arrays, and 
                dataframe for each compartment
                self.comp_arr[a].update_arrays()
                #df_sim[comp_arr[a].name] = comp_arr[d].get_df_array()"""

    def calc_tau(self):
        g_net = gna + gcl + gk + gx +g_na_k_atpase +g_kcc2
        tau = 1/g_net * cm
        return tau





    def save_to_file(self):

        with h5py.File(self.file_name, mode='a') as self.hdf:
            for i in range(len(self.comp_arr)):
                group = self.hdf.get('COMPARTMENTS')
                subgroup = group.get(self.comp_arr[i].name)
                data_array = self.comp_arr[i].get_array(self.run_t)
                subgroup.create_dataset(name=str(self.steps), data=data_array)

            for j in range(len(self.ed_arr)):
                group = self.hdf.get('ELECTRODIFFUSION')
                subgroup = group.get(self.ed_arr[j].name)
                data_array = self.ed_arr[j].ed_change_arr
                subgroup.create_dataset(name=str(self.steps), data=data_array)
