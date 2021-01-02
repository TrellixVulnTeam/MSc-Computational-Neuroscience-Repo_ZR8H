# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:55:35 2021

Calculation of the J-drift + J-diffusion for each compartment. 
Analgous to the diffusion class from Kira and Chris

@author: eshor
"""

from common import T, RTF
from compartment import Compartment
from constants import k,q,valence



class Electrodiffusion():
    """ Class to manage all the 1 dimensional electrodiffusion calculations"""
    
    def __init__ (self, comp_a: Compartment, comp_b: Compartment, d_ion:dict):
        """Creating a connection between 2 compartments
        comp_a refers to the first compartment
        comp_b refers to the second compartment
        d_ion is a dictionary in the form {'ion': Diffusion coefficient} e.g. {'na':1.33e-7}
        """
        
        self.name = comp_a.name +  ' -> ' + comp_b.name
        self.comp_a = comp_a
        self.comp_b = comp_b
        self.dx =self.comp_a.length/2 + self.comp_b.length/2
        self.ion = d_ion.keys()
        self.diff_coeff = d_ion.values()
        
        
    def calc_diffusion(self):
        """
        Calculates Fick's law for Diffusion
        F = -D*dc/dx
        """
        dc = self.comp_a.conc_i(self.ion)+ self.comp_b.conc_i(self.ion) 
        J_diffusion = -1* self.diff_coeff * dc / self.dx
        
        
    def calc_drift(self):
        """
        Calculates Ohm's law for Drift
        Drift = -D*z/RTF*dV/dx*[C]
        """
        dV = self.comp_a.V - self.comp_b.V
        z = valence(self.ion)
        J_drift = - (self.diff_coeff / RTF * z * dV / self.dx) * (self.comp_a.conc_i(self.ion)+ self.comp_b.conc_i(self.ion))
        ## Chris and Kira have a odd way of calculating the difference in concentration
        return J_drift  
    
    def calc_electrodiffusion(self):
        """Incorporates both diffusion and drift and returns an answer in Molar/s as a vector"""
        J = self.calc_drift() + self.calc_diffusion()
        return J
