# -*- coding: utf-8 -*-
"""
Created on Tue May 19 16:25:34 2020

Here I use Kira's analytical model equations to make calculating model values faster


@author: eshor
"""
import numpy as np

class Class_Single_Compartment_Analytical():
    
    R = 8.31466           # J/(Kelvin.mol) #Universal Gas Constant
    F = 96485.0           #C/mol # Faraday Constant in Volts
    T = 310.15            #Kelvin # Absolute temperature (37C)
    RTF = R*T/F           # J/C
    osmotic_perm = 0.018  #dm/s #Osmotic permeability
    partial_molar_vol_water = 0.018 #dm3/mol #Partial molar volume of water
    
    def __init__(self):
        self.vm = 0
        self.e_k = 0
        self.e_cl =0
        
    def Set_External_Concentrations(self,na=145e-3,k=3.5e-3,cl=119e-3,impermeants=29.5e-3,impermeant_charge=-0.85):
     
      self.na_out = na
      self.k_out = k
      self.cl_out = cl
      self.x_out = impermeants
      self.x_out = -1*(self.cl_out-self.na_out-self.k_out)
      self.z_out = impermeant_charge  
      self.osmol_out = self.na_out+self.k_out+self.cl_out+self.x_out
      self.nhp =0
      
      
    def Set_Leak_Conductances(self, g_na=2e-3, g_k=7e-3 , g_cl=2e-3, g_x=0):
    
      self.g_na = g_na/self.F
      self.g_k = g_k/self.F
      self.g_cl = g_cl/self.F
      self.g_x = g_x/self.F
        
         
    def Set_NaKATPase_Properties(self,pump_rate=-1):
        
       self.p_rate = pump_rate
       self.p = (10**self.p_rate)/self.F
       self.default_P=-40456
       self.jp = 10**(self.default_P/10000.0) ### Not too sure why Kira uses this Default P 
       
       
    def Set_KCC2_Properties(self,KCC2_conductance =2e-3 ):
       ''' Set the KCC2 properties (conductance in F)'''
       self.g_kcc2 = KCC2_conductance/self.F
       self.beta = self.g_k*self.g_cl - self.g_kcc2*self.g_cl +self.g_k*self.g_kcc2
       
    def Calc_Theta(self):
        
        A = -self.z_out*(self.osmol_out +self.nhp)
        B = (self.z_out**2)*(self.osmol_out+self.nhp)**2
        C = 4*(1-self.z_out**2) * self.cl_out * np.exp(-2*self.jp*self.F*self.g_kcc2 / (self.R*self.T*self.beta))
        D = self.na_out *np.exp(-3*self.jp*self.F / (self.R*self.T*self.g_na))
        D = D + self.k_out * np.exp(2*self.jp*self.F*(self.g_cl+self.g_kcc2)/(self.R*self.T*self.beta))
        numerator = A + np.sqrt(B + C*D)
        
        M = 2*(1-self.z_out) 
        N = self.na_out* np.exp(-3*self.jp*self.F/(self.R *self.T* self.g_na))
        O = self.k_out* np.exp(2*self.jp*self.F*(self.g_cl+self.g_kcc2)/(self.R*self.T*self.beta))
        denominator = M*(N+O)
        
        self.theta = numerator/denominator
        
    def Calc_Vm(self):
        
        self.vm=(-np.log(self.theta))*self.R*self.T/self.F
    
    def Calc_EK(self):
        
        self.Ek = self.vm - 2*self.jp*(self.g_cl + self.g_kcc2)/self.beta
        
    def Calc_ECl(self):
        
        self.ECl= (self.g_cl*self.vm + self.g_kcc2*self.Ek)/(self.g_cl+self.g_kcc2)
    
