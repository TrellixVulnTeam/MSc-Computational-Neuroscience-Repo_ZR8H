"""

Control the functioning of the Neural Physiological Emulator (PHANTOM)


"""

import compartment
import electrodiffusion
from common import F


class simulator:

    def __init__(self):
        """ Compartments array needs to be in the format of compartment class"""

        self.t =0 #RunTime
        self.ED_on =True
        self.constant_j_atp, self.constant_ar=False, False
        self.ed_arr = []
        self.comp_arr = []
        self.external_xflux_setup,self.xflux_setup,self.zflux_setup = True,True,True
        self.na_o, self.k_o, self.cl_o, self.x_o, self.z_o,self.osm_o = 0,0,0,0,0,0
        self.p = 0
        self.total_t, self.dt =0,0


    def add_compartment(self,comp=compartment):
        """Every compartment created needs to be added to the simulator"""
        self.comp_arr.append(comp)

    def set_electrodiffusion_properties(self,ED_on=True):
        self.ED_on = ED_on
        self.ed_arr = [electrodiffusion.Electrodiffusion(self.comp_arr[e], self.comp_arr[e + 1]) for e in range(len(self.comp_arr)-1)]

    def set_external_ion_properties(self, na_o=145e-3, k_o=3.5e-3, cl_o=119e-3, x_o=29.5e-3, z_o=-0.85):
        """
        Capacity to change the extracellular bath properties before the simulation
        """
        self.na_o, self.k_o, self.cl_o, self.x_o, self.z_o = na_o, k_o, cl_o, x_o, z_o
        self.x_o = -1 * (self.cl_o - self.na_o - self.k_o)
        self.osm_o = self.x_o + self.na_o + self.cl_o + self.k_o


    def set_j_atp(self, constant_j_atp=False,p=-1):
        self.constant_j_atp = constant_j_atp
        self.p = (10 ** p) / F


    def set_area_scale(self,constant_ar= False):
        self.constant_ar = constant_ar

    def set_timing(self, total_t,time_step):
        self.total_t, self.dt  = total_t, time_step

    def set_xflux(self, all_comps=False, comps=None,type = 'dynamic', start_t=0,end_t=0,x_conc=1e-3, flux_rate=1, z=-0.85):
        """

        @param all_comps:
        @type comps: list of compartment names
        """
        if comps is None:
            comps = []
        if all_comps:
            for a in range(len(self.comp_arr)):
                self.comp_arr[a].xflux_switch=True

        else:
            for i in range(len(comps)):
                for j in range(len(self.comp_arr)):
                    if comps[i] == self.comp_arr[j].name:
                        self.comp_arr[j].zflux_switch = True

    def set_zflux(self, all_comps=False, comps=None):
        """
        : param comps takes a list of compartment names
        """
        if comps is None:
            comps = []
        if all_comps:
            for a in range(len(self.comp_arr)):
                self.comp_arr[a].zflux_switch=True

        else:
            for i in range(len(comps)):
                for j in range(len(self.comp_arr)):
                    if comps[i]== self.comp_arr[j].name:
                        self.comp_arr[j].zflux_switch = True

    def external_xflux(self, run_t=0, start_t=0, end_t=50, xo_conc=1e-3, z=-0.85):
        """
        CHANGE THE FLUX OF IMPERMEANTS OUTSIDE THE COMPARTMENT

        """
        if start_t <= run_t <= end_t:

            if self.external_xflux_setup:
                # starting values for external x flux
                self.xo_start = self.x_o
                self.cl_o_start = self.cl_o
                self.d_xoflux = 0
                self.xo_final = self.x_o + xo_conc
                self.t_xoflux = 0
                self.xoflux_points = (end_t - start_t) * (1 / self.dt)
                self.dt_xoflux = 4 / self.xoflux_points
                self.xo_alpha = 1
                self.xo_beta = -1

            if self.x_o <= self.xo_final:
                self.external_xflux_setup = False
                self.d_xoflux = self.xo_alpha - np.e ** (self.xo_beta * self.t_xoflux)
                self.xoflux = self.d_xoflux * xo_conc
                self.x_o = self.xo_start + self.xoflux
                self.cl_o = self.cl_o_start - self.xoflux  # balancing the charges added externally
                self.t_xoflux += self.dt_xoflux

        else:
            self.xoflux = 0
            return

    #def x_flux(self):
    #def z_flux(self):
    #def xo_flux(self):



    def run_simulation(self):

        while self.run_t < self.total_t:

            if self.ED_ON:
                self.ed_dict_arr, self.ed_conc_changes_arr = [], []  # array of all the electro-diffusion dictionaries (constantly changing)

                for a in range(len(self.comp_arr)):

                    self.comp_arr[a].step(self.dt,
                                          self.na_o, self.k_o,self.cl_o,
                                          constant_j_atp=self.constant_j_atp,
                                          p = self.p)

                    )  # step for each compartment

                    if self.comp_arr[a].xflux_on:
                        self.comp_arr[a].x_flux(xflux_type, run_t, xflux_start_t, xflux_end_t, xflux_conc,
                                                       xflux_rate, xflux_charge)
                    if self.comp_arr[a].zflux_on:
                        self.comp_arr[a].z_flux(zflux_start_t, zflux_end_t, zflux_charge)

                    if self.comp_arr[a].xoflux_on:
                        self.comp_arr[a].external_xflux(run_t, xoflux_start_t, xoflux_end_t, xoflux_conc,
                                                                 xoflux_charge)
                    self.ed_dict_arr.append(self.comp_arr[a].get_ed_dict())  # electrodiffusion dictionary for each compartment

                for b in range(len(self.comp_arr) - 1):
                        self.ed_conc_changes_arr.append(self.ed_arr[b].calc_ed(self.dt, self.ed_dict_arr[b], self.ed_dict_arr[b + 1]))  # makes an array of all the ED conc changes

                for c in range(len(comp_arr) - 1):
                        comp_arr[c].ed_update(ed_conc_changes_arr[c],
                                              "positive")  # appending the electrodiffusion concentrations for each compartment
                        comp_arr[c + 1].ed_update(ed_conc_changes_arr[c], "negative")

                for d in range(len(comp_arr)):
                        comp_arr[
                            d].update_volumes()  # updates of the volumes, arrays, and dataframe for each compartment
                        if run_t != 0:
                            comp_arr[d].update_arrays()
                        df_sim[comp_arr[d].name] = comp_arr[d].get_df_array()

                if run_t != 0:
                        t_arr.append(run_t)

                    run_t += dt
                    prg.value += dt
                    lbl_prg.value = "Percent complete:" + str(round(run_t / total_t * 100, 5)) + '%'

                else:  # if you want to run with normal diffusion not ED
                    for a in range(len(comp_arr)):
                        comp_arr[a].step(dt, total_t, run_t, constant_j_atp=j_atpase_constant)
                        comp_arr[a].x_flux(xflux_start_t, xflux_end_t, xflux_conc, xflux_charge)
                        comp_arr[a].update_volumes(
                            ar_constant)  # updates of the volumes, arrays, and dataframe for each compartment
                        comp_arr[a].update_arrays()
                        df_sim[comp_arr[a].name] = comp_arr[d].get_df_array()"""

