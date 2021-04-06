# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:45:32 2021

@author: eshor

Class which defines the compartments object and related methods.

Class:
    Compartment : New compartment
Methods:
    __int__ : Initializes compartment object
    set_ion_properties: define intracellular ionic properties of the compartment (extracellular properties are imported)
    step: actions to take at every time point in the simulation
    update_volumes: change the volume of the compartment (as well as ion concentrations based on new volume)
    update_arrays: update the arrays for each parameter of the compartment
    ed_update: make changes to the compartment based on the results of electrodiffusion
    get_ed_dict: sends the current status of the compartment back to the simulation to be evaluated by electrodiffusion equations
    get_fin_vals: sends the final values of the compartment back to the simulation
    get_df_array: sends the dataframe arrays back to the simulation



"""

##################################################################################
# IMPORTS

import numpy as np

from common import oso, nao, ko, clo, xo,\
    gk, gna, gcl,gx, \
    pw, vw, RTF
from constants import F


##################################################################################
# COMPARTMENT CLASS

class Compartment():

    def __init__(self, compartment_name, radius=5e-5, length=10e-5, cm=2e-4, pkcc2=2e-3 / F, p=-1):

        self.name = compartment_name
        self.radius = radius  # in dm
        self.length = length  # in dm
        self.vw = vw
        self.pw = pw
        self.w = np.pi * (self.radius ** 2) * self.length
        self.dw = 0
        self.w2 = 0
        self.w_temp = self.w
        self.sa = 2 * np.pi * self.radius * self.length
        self.ar = self.sa / self.w
        self.C = cm
        self.FinvCAr = F / (self.C * self.ar)
        self.p_kcc2 = pkcc2
        self.j_kcc2 = 0
        self.p = (10 ** p) / F
        self.j_p = 0
        self.constant_j_p_rate = 1.1788853299370232e-09  # steady state value of model
        self.v = 0
        self.E_cl = 0
        self.E_k = 0
        self.drivingf_cl = 0

        self.na_i = 0
        self.k_i = 0
        self.cl_i = 0
        self.x_i = 0
        self.z_i = 0
        self.na_o = 0
        self.k_o = 0
        self.cl_o = 0
        self.x_o =0
        self.osm_i = 0
        self.osm_o = 0

        self.x_default = 154.962e-3

        self.z_default = -0.85

        self.na_ramp = 0
        self.diff = 0

        self.xflux_setup = True
        self.zflux_setup =True

        self.external_xflux_setup = True

        self.xflux_switch = False #if this x-flux will occur as specified
        self.zflux_switch = False


        self.xflux = 0
        self.xoflux =0

        # Zeroing Delta values
        self.d_na_i = 0
        self.d_k_i = 0
        self.d_cl_i = 0
        self.d_x_i = 0

        # Zeroing conductances
        self.g_x = 0  # basically 0 ... therefore impermeant
        self.g_na = 0
        self.g_k = 0
        self.g_cl = 0

        self.dt = 0
        self.syn_t_off = 0
        self.syn_t_on = 0

        # Zeroing arrays
        self.na_arr = []
        self.k_arr = []
        self.cl_arr = []
        self.x_arr = []

        self.z_arr = []
        self.d_na_arr = []
        self.d_k_arr = []
        self.d_cl_arr = []
        self.j_p_arr = []
        self.v_arr = []
        self.d_w_arr = []
        self.E_k_arr = []
        self.E_cl_arr = []
        self.drivingf_cl_arr = []
        self.w_arr = []
        self.ar_arr = []
        self.osm_i_arr = []
        self.osm_o_arr = []
        self.xflux_arr = []
        self.zflux_arr = []
        self.xo_arr = []

    def set_ion_properties(self, na_i=14.002e-3, k_i=122.873e-3, cl_i=5.163e-3, x_i=154.962e-3, z_i=-0.85, g_x=0e-9):
        """
        - Internal ionic concentrations and leak conductance.
        - Adjustment of starting concentrations to ensure starting electroneutrality
        - External ionic concentrations (based on imports)
        - Ionic conductances (based on imports)
        """

        self.na_i = na_i
        self.k_i = k_i
        self.cl_i = cl_i
        self.z_i = z_i
        self.x_i = x_i

        self.x_start = x_i
        self.z_start = z_i

        """if cl_i == 0:
            # setting chloride that is osmo- and electro-neutral initially.
            self.cl_i = (oso + (self.na_i + self.k_i) * (1 / self.z_i - 1)) / (1 + self.z_i)
      
           if self.k_i == 0:
            self.x_i = 155.858e-3
            self.k_i = self.cl_i-self.z_i*self.x_i-self.na_i
           else:
            self.x_i = (self.cl_i - self.k_i - self.na_i) / self.z_i """

        # Extracellular ion properties:
        self.na_o = nao
        self.k_o = ko
        self.cl_o = clo
        self.x_o = -1 * (self.cl_o - self.na_o - self.k_o)
        #self.x_o = xo
        self.osm_o = self.x_o + self.na_o + self.cl_o + self.k_o


        # Ionic conductance
        self.g_x = g_x   # basically 0 ... therefore impermeant
        self.g_na = gna
        self.g_k = gk
        self.g_cl = gcl

        #Temp values for anion fluxes:
        #self.x_ratio = 0.98  # ratio used in x_flux calculations as per Kira
        #self.x_temp_high = self.x_i(1-self.x_ratio)
        #self.x_temp_low =self.x_i(self.x_ratio)

    def set_external_ion_properties(self,na_o = 145e-3, k_o= 3.5e-3 ,cl_o = 119e-3, x_o = 29.5e-3, z_o = -0.85):
        """
        Capacity to change the extracellular bath properties before the simulation
        """

        self.na_o = na_o
        self.k_o = k_o
        self.x_o = x_o
        self.cl_o = cl_o
        self.z_o = z_o

    def osmol_neutral_start(self):
        """
        Function to ensure that the start of the simulation is osmoneutral.
        Therefore can't just start the simulation with large differences between positive and negative values
        """
        """ERAN START RAMP:
        if self.x_i * self.z_i <= -140e-3:  # -140 is a value chosen for when the ramp should be activated based on when system was crashing based on too much impermeants.
            self.diff = (self.x_i * self.z_i) + 140e-3
            self.na_ramp = self.diff * -1  # Adding some extra intracellular sodium to the compartment to offset the change of impermeants
            self.na_i += self.na_ramp"""

        """KIRA START RAMP:"""
        self.k_i = self.cl_i - self.z_i * self.x_i - self.na_i

    def step(self, dt=1e-3, total_t=120, t=0, constant_j_atp=False):
        """
        Perform a time step for the specific compartment.
        1)  Reset deltas to zero
        2)	Calculate the voltages in each compartment
        3)	Calculate the KCC2 and ATPase pump rate in each compartment (different)
        4)	Solve the ion flux in the Y direction across the membrane for each compartment (electrochemical gradient + pumps)
        5)  Update ion concentrations

        """
        # 1) Zeroing deltas
        self.dt = dt
        self.total_t = total_t
        self.t = t
        self.d_na_i = 0
        self.d_k_i = 0
        self.d_cl_i = 0
        self.d_x_i = 0

        # 2) Updating voltages
        self.v = self.FinvCAr * (self.na_i + self.k_i + (self.z_i * self.x_i) - self.cl_i)

        if self.cl_i < 0:
            print("Cl_i = " + str(self.cl_i))
            print("d_Cl_i = " + str(self.d_cl_arr[-1]))
            raise Exception("chloride log can't have a negative number")

        self.E_k = -1 * RTF * np.log(self.k_i / self.k_o)
        self.E_cl = RTF * np.log(self.cl_i / self.cl_o)
        self.drivingf_cl = self.v - self.E_cl

        # 3) Update ATPase and KCC2 pump rate
        if constant_j_atp == False:
            self.j_p = self.p * (self.na_i / nao) ** 3
        elif constant_j_atp == True:
            self.j_p = self.constant_j_p_rate

        self.j_kcc2 = self.p_kcc2 * (self.E_k - self.E_cl)

        # 4) Solve ion flux equations for t+dt from t

        self.d_na_i = - self.dt * self.ar * (
                self.g_na * (self.v + RTF * np.log(self.na_i / self.na_o)) + 3 * self.j_p)

        self.d_k_i = - self.dt * self.ar * (
                self.g_k * (self.v + RTF * np.log(self.k_i / self.k_o)) - 2 * self.j_p - self.j_kcc2)

        self.d_cl_i = + self.dt * self.ar * (self.g_cl * (self.v + RTF * np.log(self.cl_o / self.cl_i)) + self.j_kcc2)

        #self.d_x_i = - self.dt * self.ar * self.z_i * (self.g_x *
         #                                              (self.v - (RTF/self.z_i*np.log(self.x_o/self.x_temp))))

        # 5) Update ion concentrations
        self.na_i = self.na_i + self.d_na_i
        self.k_i = self.k_i + self.d_k_i
        self.cl_i = self.cl_i + self.d_cl_i

    def update_volumes(self):
        """ Calculates the new compartment volume (dm3)
        Elongation should occur radially
        """
        self.osm_i = self.na_i + self.k_i + self.cl_i + self.x_i
        self.osm_o = self.na_o + self.k_o + self.cl_o + self.x_o
        self.dw = self.dt * (self.vw * self.pw * self.sa * (self.osm_i - self.osm_o))
        self.w2 = self.w + self.dw

        self.na_i = self.na_i * self.w / self.w2
        self.k_i = self.k_i * self.w / self.w2
        self.cl_i = self.cl_i * self.w / self.w2
        self.x_i = self.x_i * self.w / self.w2

        self.w = self.w2

        self.ar = self.sa / self.w
        self.FinvCAr = F / (self.C * self.ar)

    def update_arrays(self):
        """
        Update arrays such that they reflect the actual values
        """
        self.na_arr.append(self.na_i * 1000)
        self.k_arr.append(self.k_i * 1000)
        self.cl_arr.append(self.cl_i * 1000)
        self.x_arr.append(self.x_i * 1000)
        self.z_arr.append(self.z_i)
        self.w_arr.append(self.w * (10 ** 12))
        self.ar_arr.append(self.ar)
        self.v_arr.append(self.v * 1000)
        self.d_na_arr.append(self.d_na_i * 1000)
        self.d_k_arr.append(self.d_k_i * 1000)
        self.d_cl_arr.append(self.d_cl_i * 1000)
        self.j_p_arr.append(self.j_p)
        self.d_w_arr.append(self.dw * 1000)
        self.E_k_arr.append(self.E_k * 1000)
        self.E_cl_arr.append(self.E_cl * 1000)
        self.drivingf_cl_arr.append(self.drivingf_cl * 1000)
        self.osm_i_arr.append(self.osm_i * 1000)
        self.osm_o_arr.append(self.osm_o * 1000)
        self.xflux_arr.append(self.xflux * 1000)
        self.xo_arr.append(self.x_o*1000)

    def ed_update(self, ed_change: dict, sign="positive"):
        """
        Receives a dictionary and update
        """
        if sign == "positive":
            self.na_i += ed_change["na"]
            self.cl_i += ed_change["cl"]
            self.k_i += ed_change["k"]
            self.x_i += ed_change["x"]
        elif sign == "negative":
            self.na_i -= ed_change["na"]
            self.cl_i -= ed_change["cl"]
            self.k_i -= ed_change["k"]
            self.x_i -= ed_change["x"]

    def get_ed_dict(self):
        ed_dict = {"na": self.na_i, "k": self.k_i, "cl": self.cl_i, "x": self.x_i, "Vm": self.v}
        return ed_dict

    def get_fin_vals(self):

        valstring = self.name + " final values: "
        valstring = valstring + " Na: " + str(self.na_arr[-1])
        valstring = valstring + " K: " + str(self.k_arr[-1])
        valstring = valstring + " Cl " + str(self.cl_arr[-1])
        valstring = valstring + " X " + str(self.x_arr[-1])
        valstring = valstring + " Vm: " + str(self.v_arr[-1])
        valstring = valstring + " Volume: " + str(self.w_arr[-1])

        return valstring

    def get_df_array(self):
        df_arr = [self.radius, self.length, self.w, self.na_i, self.k_i, self.cl_i, self.x_i, self.z_i, self.p,
                  self.p_kcc2, self.v, self.E_k, self.E_cl, ]
        return df_arr

    def add_synapse(self, syn_t_on=0, syn_t_off=0, Inhibitory_synapses=0, Excitatory_synapses=0, Ramp=True):
        """
        Set time synapse is turned on and off
        Inhibitory_synapses = number of inhibitory synapses onto compartment.
        Inhibition results in  Chloride conductance change (Simulating GABA-R activation)

        Excitatory_synapses = number of excitatory synapses onto compartment.
        Excitation results in  Sodium conductance change (Simulating NMDA-R activation)

        """
        self.syn_t_on = syn_t_on
        self.syn_t_off = syn_t_off

        if Inhibitory_synapses != 0:
            self.g_cl = self.g_cl * 10 * Inhibitory_synapses
        if Excitatory_synapses != 0:
            self.g_na = self.g_na * 10 * Excitatory_synapses

        return

    def x_flux(self, run_t=0, start_t=0, end_t=50, x_conc=1e-3, z=-0.85):
        """
        FLUX IMPERMEANTS INTO THE COMPARTMENT
        :param start_t:  Start time to impermeant flux
        :param end_t:  End time of impermeant fulx
        :param x_conc: Concentration of impermenat to add to the model
        :param z: charge of impermeant to add to the model
        :return:
        """

        if (start_t <= run_t <= end_t) and self.xflux_switch:

            if self.xflux_setup:

                #starting values for flux
                self.x_start = self.x_i
                self.d_xflux = 0
                self.x_final = self.x_i + x_conc
                self.t_xflux = 0
                self.flux_points = (end_t - start_t) * (1/self.dt)
                self.dt_xflux = 4/self.flux_points
                self.alpha = 1
                self.beta = -1

            if self.x_i <=self.x_final:
                self.xflux_setup = False
                self.d_xflux = self.alpha - np.e**(self.beta * self.t_xflux)
                self.xflux = self.d_xflux * x_conc
                self.x_i = self.x_start+self.xflux
                self.t_xflux += self.dt_xflux

        else:
            self.xflux = 0
            return

    def z_flux(self, start_t=0, end_t=50, z=-0.85):
        """

        :param start_t:  Start time to impermeant flux
        :param end_t:  End time of impermeant fulx
        :param z: charge of impermeant to add to the model
        :return:
        """
        if self.zflux_setup:
            self.z_diff = z-self.z_i
            #self.x_mol_start =  self.w * self.x_i


        if (start_t <= self.t <= end_t) and self.zflux_switch:

            t_diff = (end_t - start_t) / self.dt

            z_inc = self.z_diff / t_diff
            self.z_i += z_inc
            #self.x_i = self.x_mol_start/self.w
            self.zflux_setup = False
        else:
            z_inc = 0
            return

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
                self.cl_o = self.cl_o_start - self.xoflux #balancing the charges added externally
                self.t_xoflux += self.dt_xoflux

        else:
            self.xoflux = 0
            return

    def get_x_value(self, x_i_conc=154.962e-3, t_current=0, t_total=0):

        t_ramp_on = t_total / 10
        t_ramp_off = t_total / 5 * 3

        if t_current <= t_ramp_on:
            x_value = self.x_default
            return x_value

        elif t_current >= t_ramp_off:
            x_value = x_i_conc
            return x_value

        else:
            time_diff = t_ramp_off - t_ramp_on
            conc_diff = self.x_start - self.x_default
            grad = conc_diff / time_diff
            t_ramp = t_current - t_ramp_on
            x_value = self.x_default + t_ramp * grad
            return x_value
