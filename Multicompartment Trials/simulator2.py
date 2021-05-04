"""

Control the functioning of the Neural Physiological Emulator (PHANTOM)


"""

import compartment
import electrodiffusion
from common import F
import numpy as np
import time
import pandas as pd
import h5py


class simulator:

    def __init__(self, file_name =""):
        """ Compartments array needs to be in the format of compartment class"""
        self.file_name = file_name
        self.file_name = "\ "+ file_name

        try:
            with h5py.File(self.file_name, mode='w') as self.hdf:
                self.hdf.create_group('COMPARTMENTS')
                self.hdf.create_group('ELECTRODIFFUSION')
                self.hdf.create_group("TIMING")
                self.hdf.create_group("X-FLUX")
                self.hdf.create_group("SYNAPSES")
                print("simulation file ('"+file_name+"') created")


        except:
            raise ("file not created")

        # COMPS = hdf.put('COMPARTMENTS/COMP1', df_comp1, format='table')
        # ELECTRODIFFUSION = hdf.create_group('ELECTRODIFFUSION')




        self.num_comps = 0
        self.run_t = 0
        self.start_t = time.time()
        self.one_percent_t = 0
        self.end_t = 0
        self.interval_num = 1


        self.output_arr = []
        self.output_intervals = []
        self.ED_on = True
        self.ed_dict_arr, self.ed_conc_changes_arr = [], []
        self.constant_j_atp, self.constant_ar = False, False
        self.ed_arr = []
        self.comp_arr = []
        self.external_xflux_setup, self.xflux_setup, self.zflux_setup = True, True, True
        self.na_o, self.k_o, self.cl_o, self.x_o, self.z_o, self.osm_o = 0, 0, 0, 0, 0, 0
        self.p = 0
        self.total_t, self.dt = 0, 0
        self.xoflux_switch = False
        self.xoflux_params = {"start_t": 0, "end_t": 0, "xo_conc": 0, "zo": 0}
        self.xoflux_setup = True
        self.xo_start, self.cl_o_start, self.d_xoflux, self.xo_final, self.xo_flux, self.t_xoflux = 0, 0, 0, 0, 0, 0
        self.xoflux_points, self.dt_xoflux, self.xo_alpha, self.xo_beta = 0, 0, 0, 0

    def add_compartment(self, comp=compartment):
        """Every compartment created needs to be added to the simulator"""

        new_comp = comp.get_array(time=0)
        with h5py.File(self.file_name, mode='a') as self.hdf:
            group = self.hdf.get('COMPARTMENTS')
            subgroup = group.create_group(name=comp.name)
            subgroup.create_dataset(name='0', data=new_comp)

        self.num_comps += 1
        self.comp_arr.append(comp)

    def set_electrodiffusion_properties(self, ED_on=True):
        self.ED_on = ED_on
        self.ed_setup_arr = []
        with h5py.File(self.file_name, mode='a') as self.hdf:
            comp_group = self.hdf.get('COMPARTMENTS')
            ed_group = self.hdf.get('ELECTRODIFFUSION')

            for e in range(self.num_comps - 1):
                name = self.comp_arr[e].name + ' <-> ' + self.comp_arr[e + 1].name
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
                self.ed_setup_arr.append(ed.ed_setup)

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
            timing.create_dataset("DT", data=self.dt)
            timing.create_dataset("TOTAL_T", data=self.total_t)

        self.total_steps = self.total_t / self.dt

        self.output_intervals = [0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75, 1]
        self.output_arr = [self.output_intervals[a] * self.total_t for a in range(len(self.output_intervals))]


        self.interval_step =self.total_t/intervals
        self.interval_arr = [self.interval_step*i for i in range(intervals)]
       # print(self.interval_arr)

    def set_xflux(self, all_comps=False, comps=None, type='dynamic', start_t=0, end_t=0, x_conc=1e-3, flux_rate=1,
                  z=-0.85):

        xflux_params = []
        xflux_comps = []

        for i in range(len(comps)):
            for j in range(len(self.comp_arr)):
                if comps[i] == self.comp_arr[j] or all_comps:
                    xflux_comps.append(self.comp_arr[j])
                    params = []  # array of the xflux settings for a compartments
                    params.append(j)  # compartment number
                    if type == "dynamic":
                        params.append(0)  # whether xflux is static or dynamic
                    elif type == "static":
                        params.append(1)
                    params.append(start_t)
                    params.append(end_t)
                    params.append(x_conc)
                    params.append(z)
                    params.append(flux_rate)
                    xflux_params.append(params)

        with h5py.File(self.file_name, mode='a') as self.hdf:
            xflux_group = self.hdf.get("X-FLUX")
            xflux_setup = xflux_group.create_dataset("X-FLUX SETUP",data=xflux_params)
            for c in range(len(xflux_comps)):
                xflux_group.create_dataset("X-FLUX " + xflux_comps[c],data=[])



    # xflux_params is a dictionary sent to compartments that have the xflux_switch on

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


    def run_simulation(self):

        while self.run_t < self.total_t:

            if self.ED_on:
                self.ed_dict_arr =[]
                self.ed_conc_changes_arr = []

                for a in self.gen_comps(self.comp_arr):

                    a.step(dt= self.dt,
                           na_o= self.na_o, k_o= self.k_o, cl_o= self.cl_o,
                           constant_j_atp=self.constant_j_atp,
                           p=self.p, p_kcc2 = 2e-3 / F)


                    # step for each compartment

                    if a.xflux_switch and \
                            (a.xflux_params["start_t"] <= self.run_t <= a.xflux_params[
                                "end_t"]):
                        a.x_flux()

                    if a.zflux_switch and \
                            (a.zflux_params["start_t"] <= self.run_t <= a.zflux_params[
                                "end_t"]):
                        a.z_flux()

                    if self.xoflux_switch and \
                            self.xoflux_params["start_t"] <= self.run_t <= self.xoflux_params["end_t"]:
                        self.xoflux()

                    self.ed_dict_arr.append(a.get_ed_dict())
                    # electrodiffusion dictionary for each compartment

                for b in range(len(self.ed_arr)):
                    self.ed_conc_changes_arr.append(self.ed_arr[b].calc_ed(self.dt, self.ed_dict_arr[b],
                                                                           self.ed_dict_arr[b + 1]))

                    # makes an array of all the ED conc changes
                self.c2 = 0
                for c in self.gen_comps(self.comp_arr[0:-2]):
                    c.ed_update(self.ed_conc_changes_arr[self.c2],
                                "positive")
                    c.ed_update(self.ed_conc_changes_arr[self.c2], "negative")
                    self.c2 = + 1
                    # appending the electrodiffusion concentrations for each compartment

                for d in self.gen_comps(self.comp_arr):
                    d.update_volumes(self.dt, self.osm_o,
                                     self.constant_ar)  # updates of the volumes, arrays, and dataframe for each compartment




                for f in range(len(self.output_arr)):
                    if round(self.run_t, 7) == self.output_arr[f]:
                        if f == 2:
                            self.one_percent_t = time.time() - self.start_t
                            self.hundred_percent_t = self.one_percent_t * 100
                            print(str(self.output_intervals[f] * 100) + " % complete in " + str(
                                round(self.one_percent_t, 2)) + " s")
                            print("Estimated time to complete :" + str(round(self.hundred_percent_t / 60, 2)) + " minutes")
                        else:
                            print(str(self.output_intervals[f] * 100) + " % complete in " + str(round(time.time() - self.start_t, 2)) + " s")

            self.run_t += self.dt
            if round(self.run_t,6) == round(self.interval_num * self.interval_step,6):
                self.interval_num += 1
                self.save_to_file(self.run_t)
        # print(str(round(self.run_t/self.total_t*100,2))+"%")

        self.end_t = time.time()


    """elif not self.ED_on:  # if you want to run with normal diffusion not ED
            for a in range(len(self.comp_arr)):
                self.comp_arr[a].step()
                self.comp_arr[a].x_flux()
                self.comp_arr[a].update_volumes(ar_constant)  # updates of the volumes, arrays, and 
                dataframe for each compartment
                self.comp_arr[a].update_arrays()
                #df_sim[comp_arr[a].name] = comp_arr[d].get_df_array()"""

    def save_to_file(self,time):

        with h5py.File(self.file_name, mode='a') as self.hdf:
            for i in range(len(self.comp_arr)):
                group = self.hdf.get('COMPARTMENTS')
                subgroup = group.get(self.comp_arr[i].name)
                data_array= self.comp_arr[i].get_array()
                subgroup.create_dataset(name=str(time), data=data_array)

            for j in range(len(self.ed_arr)):
                group = self.hdf.get('ELECTRODIFFUSION')
                subgroup = group.get(self.ed_arr[j].name)
                data_array = self.ed_arr[j].ed_change_arr
                subgroup.create_dataset(name=str(time), data=data_array)




