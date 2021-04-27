"""

Control the functioning of the Neural Physiological Emulator (PHANTOM)


"""

import compartment
import electrodiffusion
from common import F
import numpy as np
from os import system, name


class simulator:

    def __init__(self):
        """ Compartments array needs to be in the format of compartment class"""

        self.run_t = 0
        self.t_arr = []
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
        self.comp_arr.append(comp)

    def set_electrodiffusion_properties(self, ED_on=True):
        self.ED_on = ED_on
        self.ed_arr = [electrodiffusion.Electrodiffusion(self.comp_arr[e], self.comp_arr[e + 1]) for e in
                       range(len(self.comp_arr) - 1)]

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

    def set_timing(self, total_t, time_step):
        self.total_t, self.dt = total_t, time_step
        for i in range(len(self.comp_arr)):
            self.comp_arr[i].dt=time_step

    def set_xflux(self, all_comps=False, comps=None, type='dynamic', start_t=0, end_t=0, x_conc=1e-3, flux_rate=1,
                  z=-0.85):
        for i in range(len(comps)):
            for j in range(len(self.comp_arr)):
                if comps[i] == self.comp_arr[j].name or all_comps:
                    self.comp_arr[j].xflux_switch = True
                    self.comp_arr[j].xflux_params = {"type": type, "start_t": start_t, "end_t": end_t, "x_conc": x_conc,
                                                     "flux_rate": flux_rate, "z": z}

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

                for a in range(len(self.comp_arr)):

                    self.comp_arr[a].step(self.dt,
                                          self.na_o, self.k_o, self.cl_o,
                                          constant_j_atp=self.constant_j_atp,
                                          p=self.p)

                    # step for each compartment

                    if self.comp_arr[a].xflux_switch and \
                            (self.comp_arr[a].xflux_params["start_t"] <= self.run_t <= self.comp_arr[a].xflux_params[
                                "end_t"]):
                        self.comp_arr[a].x_flux()

                    if self.comp_arr[a].zflux_switch and \
                            (self.comp_arr[a].zflux_params["start_t"] <= self.run_t <= self.comp_arr[a].zflux_params[
                                "end_t"]):
                        self.comp_arr[a].z_flux()

                    if self.xoflux_switch and \
                            self.xoflux_params["start_t"] <= self.run_t <= self.xoflux_params["end_t"]:
                        self.xoflux()

                    self.ed_dict_arr.append(self.comp_arr[a].get_ed_dict())
                    # electrodiffusion dictionary for each compartment

                for b in range(len(self.comp_arr) - 1):
                    self.ed_conc_changes_arr.append(self.ed_arr[b].calc_ed(self.dt, self.ed_dict_arr[b],
                                                                           self.ed_dict_arr[b + 1]))
                    # makes an array of all the ED conc changes

                for c in range(len(self.comp_arr) - 1):
                    self.comp_arr[c].ed_update(self.ed_conc_changes_arr[c],
                                               "positive")
                    self.comp_arr[c + 1].ed_update(self.ed_conc_changes_arr[c], "negative")
                    # appending the electrodiffusion concentrations for each compartment

                for d in range(len(self.comp_arr)):
                    self.comp_arr[
                        d].update_volumes(self.dt, self.osm_o, self.constant_ar)  # updates of the volumes, arrays, and dataframe for each compartment
                    if self.run_t != 0:
                        self.comp_arr[d].update_arrays()
                    # df_sim[comp_arr[d].name] = comp_arr[d].get_df_array()


                self.t_arr.append(self.run_t)
                self.run_t += self.dt
                #print(str(round(self.run_t/self.total_t*100,2))+"%")

            """elif not self.ED_on:  # if you want to run with normal diffusion not ED
                for a in range(len(self.comp_arr)):
                    self.comp_arr[a].step()
                    self.comp_arr[a].x_flux()
                    self.comp_arr[a].update_volumes(ar_constant)  # updates of the volumes, arrays, and 
                    dataframe for each compartment
                    self.comp_arr[a].update_arrays()
                    #df_sim[comp_arr[a].name] = comp_arr[d].get_df_array()"""
