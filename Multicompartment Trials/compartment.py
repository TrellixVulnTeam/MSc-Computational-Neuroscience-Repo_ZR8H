# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:45:32 2021

@author: eshor



"""
import numpy as np



# import copy -- WTF is that
from common import default_radius, default_length, default_p, default_Cm, \
    oso, nao, ko, clo, \
    gk, gna, gcl, \
    pw, vw, RTF
from constants import F


class Compartment():

    def __init__(self, compartment_name, radius=default_radius, length=default_length, Cm=default_Cm, pkcc2=2e-3 / F,
                 p=default_p):
        self.name = compartment_name
        self.radius = radius  # in dm
        self.length = length  # in dm
        self.w = (np.pi) * (self.radius ** 2) * self.length
        self.w_temp = self.w
        self.sa = 2 * (np.pi) * (self.radius) * (self.length)
        self.ar = self.sa / self.w
        self.C = Cm
        self.FinvCAr = F / (self.C * self.ar)
        self.p_kcc2 = pkcc2
        self.p = p
        self.V=0
        self.E_cl=0
        self.E_k=0



    def set_ion_properties(self, na_i=14.002e-3, k_i=122.873e-3, cl_i=5.163e-3, x_i = 154.962e-3, z_i=-0.85, g_x=0e-9):

        # A) Intracellular ion properties:
        self.na_i = na_i
        self.k_i = k_i
        self.cl_i = cl_i
        self.z_i = z_i
        self.x_i = x_i

        """if cl_i == 0:
            # setting chloride that is osmo- and electro-neutral initially.
            self.cl_i = (oso + (self.na_i + self.k_i) * (1 / self.z_i - 1)) / (1 + self.z_i)
      
      if self.k_i == 0:
            self.x_i = 155.858e-3
            self.k_i = self.cl_i-self.z_i*self.x_i-self.na_i
      else:
            self.x_i = (self.cl_i - self.k_i - self.na_i) / self.z_i """

        self.g_x = g_x  # basically 0 ... therefore impermeant
        self.g_na = gna
        self.g_k = gk
        self.g_cl = gcl



        # C) Extracellular ion properties:
        self.na_o = nao
        self.k_o = ko
        self.cl_o = clo
        self.osm_o = oso

        # D) Zeroing Delta values
        self.d_na_i = 0
        self.d_k_i = 0
        self.d_cl_i = 0
        self.d_x_i = 0


        # F) Zeroing arrays
        self.na_arr = []
        self.k_arr = []
        self.cl_arr = []
        self.x_arr = []
        self.d_na_arr = []
        self.d_k_arr = []
        self.d_cl_arr = []
        self.v_arr = []
        self.d_w_arr = []
        self.E_k_arr = []
        self.E_cl_arr = []
        self.w_arr = []
        self.ar_arr = []
        self.osm_i_arr=[]
        self.osm_o_arr =[]

    def step(self, dt=1e-3):
        """
        Perform a time step for the specific compartment.
        
        1)	Calculate the voltages in each compartment
        2)	Calculate the KCC2 and ATPase pump rate in each compartment (different)
        3)	Solve the ion flux in the Y direction across the membrane for each compartment (electrochemical gradient + pumps)
        4)	Solve the ion flux for each ion in the x direction (electrodiffusion equation)
        5)	Make adjustments for the ion concentrations based on movement between compartments
        6)	Calculate the volume change in each compartment
        7)	Correct ion concentrations based on new volumes

        
        """

        self.dt = dt

        # 1) Updating voltages
        self.v = self.FinvCAr * (self.na_i + self.k_i + (self.z_i * self.x_i) - self.cl_i)
        self.E_k = RTF * np.log(ko / self.k_i)
        self.E_cl = RTF * np.log(self.cl_i / clo)

        # 2) Update ATPase and KCC2 pump rate
        self.j_p = self.p * (self.na_i / nao) ** 3
        self.j_kcc2 = self.p_kcc2 * (self.E_k - self.E_cl)

        # 3) Solve ion flux equations for t+dt from t

        self.d_na_i = - self.dt * self.ar * (self.g_na * (self.v - RTF * np.log(self.na_o / self.na_i)) + 3 * self.j_p)
        self.d_k_i = - self.dt * self.ar * (self.g_k * (self.v - RTF * np.log(self.k_o / self.k_i)) - 2 * self.j_p - self.j_kcc2)
        self.d_cl_i = + self.dt * self.ar * (self.g_cl * (self.v + RTF * np.log(self.cl_o / self.cl_i)) + self.j_kcc2)
        self.na_i = self.na_i + self.d_na_i
        self.k_i = self.k_i + self.d_k_i
        self.cl_i = self.cl_i + self.d_cl_i

        # 4) Electrodiffusion calculations

        # 6) Test Arrays:



    def update_volumes(self):
        # 6) Update volume
        ''' Calculates the new compartment volume (dm3)
        Elongation should occur length ways not radially
        '''
        self.osm_i = self.na_i + self.k_i + self.cl_i + self.x_i
        self.dw = self.dt * pw * vw * self.sa * (self.osm_i - self.osm_o)
        self.w2 = self.w + self.dw

        self.na_i = self.na_i * self.w / self.w2
        self.k_i = self.k_i * self.w / self.w2
        self.cl_i = self.cl_i * self.w / self.w2
        self.x_i = self.x_i * self.w / self.w2

        self.w = self.w2

        #self.length = self.w / (np.pi * self.radius ** 2)
        self.radius = np.sqrt(self.w /(np.pi*self.length))
        self.sa = 2 * (np.pi) * (self.radius) * (self.length)

        self.ar = self.sa / self.w
        self.FinvCAr = F / (self.C * self.ar)

        ##### Radius and length will also have to be updated for electrodiffusion calculation

    def update_arrays(self):

        self.na_arr.append(self.na_i*1000)
        self.k_arr.append(self.k_i*1000)
        self.cl_arr.append(self.cl_i*1000)
        self.x_arr.append(self.x_i*1000)
        self.w_arr.append(self.w*(10**12))
        self.ar_arr.append(self.ar)
        self.v_arr.append(self.v*1000)
        self.d_na_arr.append(self.d_na_i)
        self.d_k_arr.append(self.d_k_i)
        self.d_cl_arr.append(self.d_cl_i)
        self.d_w_arr.append(self.dw)
        self.E_k_arr.append(self.E_k)
        self.E_cl_arr.append(self.E_cl)
        self.osm_i_arr.append(self.osm_i)
        self.osm_o_arr.append(self.osm_o)

    def ed_update(self, ed_change: dict):

        self.na_i += ed_change["na"]
        self.cl_i += ed_change["cl"]
        self.k_i += ed_change["k"]
        self.x_i += ed_change["x"]

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
        df_arr = [self.radius,self.length,self.w,self.na_i,self.k_i,self.cl_i,self.x_i,self.z_i,self.V,self.E_k,self.E_cl,self.p,self.p_kcc2]
        return df_arr



