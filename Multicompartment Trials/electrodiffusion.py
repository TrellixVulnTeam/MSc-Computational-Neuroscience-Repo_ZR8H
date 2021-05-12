# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:55:35 2021

Calculation of the J-drift + J-diffusion for each compartment. 
Analgous to the diffusion class from Kira and Chris

@author: eshor
"""

from common import T, RTF,val,diff_constants
from compartment import Compartment
from constants import k,q,valence
import pandas as pd
import h5py


class Electrodiffusion:
    """ Class to manage all the 1 dimensional electrodiffusion calculations"""
    
    def __init__ (self, comp_a_name= '', comp_a_length=10, comp_b_name ='', comp_b_length=10):
        """Creating a connection between 2 compartments
        comp_a refers to the first compartment
        comp_b refers to the second compartment
        d_ion is a dictionary in the form {'ion': Diffusion coefficient} e.g. {'na':1.33e-7}
        """
        
        self.name = comp_a_name +  ' <-> ' + comp_b_name
        self.comp_a = comp_a_name
        self.comp_b = comp_b_name
        self.dx =comp_a_length/2 + comp_b_length/2
        #self.bound_na_arr, self.bound_k_arr,self.bound_cl_arr,self.bound_x_arr = [],[],[],[] #arr of the various changes occuring at the boundary

        self.ed_setup = [self.name,self.comp_a,self.comp_b,self.dx]



    def calc_diffusion(self,ion="",conc_a=0,conc_b=0):
        """
        Calculates Fick's law for Diffusion
        F = -D*dc/dx
        """
        d = diff_constants[ion]
        dc = conc_b - conc_a
        j_diffusion = -1 * d * dc / self.dx
        return j_diffusion
        
        
    def calc_drift(self,ion = "",conc_a=0, conc_b=0,dV=0,):
        """
        Calculates Ohm's law for Drift
        Drift = -D*z/RTF*dV/dx*[C]
        """
        z = val[ion]
        d = diff_constants[ion]
        j_drift = - (d / RTF * z * dV / self.dx) * (conc_a + conc_b)
        # Chris and Kira have a odd way of calculating the difference in concentration
        return j_drift
    
    def calc_ed(self,dt=1e-3, comp_a_ed_dict={"na":0,"k":0, "cl":0,"x":0,"Vm":0},comp_b_ed_dict={"na":0,"k":0, "cl":0,"x":0,"Vm":0}):
        """Incorporates both diffusion and drift and returns an answer in Molar/s as a vector

        * Note that the flux between compartment a to b = flux from b to a
        In Chris' electrodiffusion updates he divides the drift by 2 ... Not too sure why?

        """

        self.ed_change_dict = {"na": 0, "k": 0, "cl": 0, "x": 0}
        self.ed_change_arr = [] #array which will be saved as the dataset
        dv = comp_a_ed_dict["Vm"] - comp_b_ed_dict["Vm"]
        ions = list(self.ed_change_dict)

        for i in range(4):
            ion = ions[i]
            self.ed_change_dict[ion] += self.calc_drift(ion, comp_a_ed_dict[ion], comp_b_ed_dict[ion], dv)/2
            self.ed_change_dict[ion] += self.calc_diffusion(ion, comp_a_ed_dict[ion], comp_b_ed_dict[ion])
            self.ed_change_dict[ion] *= dt
            self.ed_change_arr.append(self.ed_change_dict[ion])

        return self.ed_change_dict


