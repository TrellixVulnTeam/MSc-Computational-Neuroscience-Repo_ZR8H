# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 16:05:53 2020

Defined here is a class which runs the single compartment neuronal model based on the pump-leak mechanism.

Method 1: Initializes the SingleCompartment_PumpLeak class.
    - 

@author: eshor
"""

import numpy as np


class SingleCompartment_PumpLeak():

  ### CONSTANTS  
  R = 8.31466           # J/(Kelvin.mol) #Universal Gas Constant
  F = 96485.0           #C/mol # Faraday Constant in Volts
  T = 310.15            #Kelvin # Absolute temperature (37C)
  RTF = R*T/F           # J/C
  osmotic_perm = 0.018  #dm/s #Osmotic permeability
  partial_molar_vol_water = 0.018 #dm3/mol #Partial molar volume of water
  default_p = 0.1 #C/dm2.s #default pump rate

  def __init__(self, radius = 5e-5, length= 25e-5, capacitance = 2e-4):
      ''' Initializes the SingleCompartment_PumpLeak class radius and length in dm; capacitance in F'''
      self.radius = radius
      self.length = length
      self.cm = capacitance
      self.volume = np.pi*self.radius**2*self.length
      self.area = 2*np.pi*self.radius*self.length
      self.area_scale = self.area/self.volume
      self.f_inv_c_ar = self.F / (self.cm *self.area_scale)


  def Set_Concentrations(self, Na = [14.002e-3, 145e-3], K = [122.873e-3, 3.5e-3], Cl = [5.163e-3, 119e-3], Impermeants = [154.962e-3, 29.5e-3], Impermeant_charges = [-0.85,-0.85] ):
      ''' Sets the internal and external ionic concentrations (Molar) '''
      self.na_in = Na[0]
      self.na_out = Na[1]
      self.k_in = K[0]
      self.k_out = K[1]
      self.cl_in = Cl[0]
      self.cl_out = Cl[1]
      self.x_in = Impermeants[0]
      self.x_out = Impermeants[1]
      self.x_out = -1*0.2*(self.cl_out-self.na_out-self.k_out)
      self.z_in = Impermeant_charges[0]
      self.z_out = Impermeant_charges[1]
      # Kira uses these to get starting extravellular impearmeants xe1=-1*(cle-nae-ke)
      # xe=xe1*0.2


  def Set_LeakConductances(self, gNa= 2e-3 , gK= 7e-3 ,gCl=2e-3 , gImpermeants=0):
      ''' Sets the conductances (Siemens) of the leak channels for each ion species '''
      self.g_na = gNa/self.F
      self.g_k = gK/self.F
      self.g_cl = gCl/self.F
      self.g_x = gImpermeants/self.F

  def Initialize_Arrays(self):
      ''' Initializes the SingleCompartment_PumpLeak class'''
      self.na_arr = []
      self.k_arr = []
      self.cl_arr = []
      self.x_arr = []
      self.vm_arr = []
      self.t_arr = []           
      self.volume_arr = []
      self.e_cl_arr =[]
      self.e_k_arr = []


  def Calc_MembraneVoltage(self):
       ''' Calculate the membrane voltage (V) based on the Charge Difference approach Fraser and Huang (2004)'''
       self.vm = self.f_inv_c_ar*(self.na_in+self.k_in-self.cl_in+(self.x_in*self.z_in))
       
  def Set_NaKATPase_Properties(self,time_off=300,time_on=900,ATP_pump_rate =-1, switch_at_start = 'on' ):
       ''' Set the Na-K ATPase properties'''
       if (switch_at_start == "on") or (switch_at_start == 1) or (switch_at_start == "ON") or (switch_at_start=="On") :
           self.atp_switch = 1
       else: 
           self.atp_switch = 0
           
       self.atp_t_on = time_on
       self.atp_t_off = time_off
       self.p_default = ATP_pump_rate
       self.p_effective = (10**self.p_default)/self.F
       
 
  def Set_KCC2_Properties(self,time_on=0,time_off=1000000,KCC2_conductance =2e-3 ):
       ''' Set the KCC2 properties (conductance in F)'''
       self.kcc2_t_on = time_on
       self.kcc2_t_off = time_off
       self.g_kcc2 = KCC2_conductance/self.F

  def Set_Timing(self,duration = 1000, time_step = 1e-3, num_intervals = 180000):
       ''' Set the timing for the simulation (seconds) '''
       self.t = 0                                       # real time in seconds
       self.total_t = duration                          # total duration
       self.dt = time_step                              # time step
       self.total_steps = round(self.total_t/self.dt)   #total time steps
       self.total_bins = num_intervals                  # number of plot points
       self.bin = 1                                     # counter for plotting 
       self.ts = self.total_t/self.total_bins           #time in between plot points   
       
  def Check_ATP_switch_position(self):
      ''' Assess if ATPase is on or off'''
      if (self.t>self.atp_t_off) & (self.t<self.atp_t_on):
          self.atp_switch = 0
      else:
          self.atp_switch = 1    
          
  def Calc_ATPase_rate(self):
      
       p_rate = self.p_default
       self.p_max = -5 
       self.ramp_change = (p_rate-self.p_max)/(12*10**4)/8   #equivalent to em in Kira's code
          
       if self.atp_switch  ==1:
              if p_rate > self.p_max:
                  p_rate -= self.ramp_change
                  self.p_effective = (10**p_rate)/self.F
       else:
              if p_rate < self.p_max:
                  p_rate += self.ramp_change
                  self.p_effective = (10**p_rate)/self.F   
                  
       self.j_atp = self.p_effective*(self.na_in/self.na_out)**3
          
  
  def Calc_new_volume(self):
      ''' Calculates the new compartment volume (dm3)'''
      self.osmol_in = self.na_in+self.k_in+self.cl_in+self.x_in  
      self.osmol_out = self.na_out + self.k_out +self.cl_out+self.x_out 
      self.radius = np.sqrt(self.volume/(np.pi*self.length))
      self.area = 2*np.pi*self.radius*self.length 
      self.dw=self.dt*self.partial_molar_vol_water*self.osmotic_perm*self.area*(self.osmol_in-self.osmol_out)  
      self.volume2 = self.volume+self.dw
          
  def Calc_new_ion_conc(self):  
      '''Update Ion concentrations based on the differential equations'''
      d_na = -self.dt*self.area_scale*(self.g_na*(self.vm-self.RTF*np.log(self.na_out/self.na_in)) + self.atp_switch*3*self.j_atp) # - (1/self.volume)*dw*self.na_in
      d_k = -self.dt*self.area_scale*(self.g_k*(self.vm-self.RTF*np.log(self.k_out/self.k_in)) - self.atp_switch*2*self.j_atp - self.j_kcc2)#  - (1/self.volume)*dw*self.k_in
      d_cl = self.dt*self.area_scale*(self.g_cl*(self.vm-self.RTF*np.log(self.cl_in/self.cl_out)) - self.j_kcc2) # - (1/self.volume)*dw*self.cl_in
      self.na_in += d_na
      self.k_in += d_k
      self.cl_in += d_cl
      
  def Adjust_concentrations_for_volume(self):
      ''' Recalculates concentrations based on updated volumes'''
      self.na_in *= self.volume/self.volume2
      self.k_in *= self.volume/self.volume2
      self.cl_in *= self.volume/self.volume2
      self.x_in *= self.volume/self.volume2
            
  def Update_arrays(self):
      ''' Update Arrays'''
      self.vm_arr.append(self.vm*1e3)
      self.e_cl_arr.append(self.e_cl*1e3)
      self.e_k_arr.append(self.e_k*1e3)
      self.na_arr.append(self.na_in*1e3)
      self.k_arr.append(self.k_in*1e3)
      self.cl_arr.append(self.cl_in*1e3)
      self.x_arr.append(self.x_in*1e3)
      self.volume_arr.append(self.volume)
      self.t_arr.append(self.t)
        
  def Simulate_PLM(self):
      ''' Simulate the pump leak model mechanism''' 
      
      for i in range(1, self.total_steps):  
          
          self.Calc_MembraneVoltage()
          self.e_k = self.RTF*np.log(self.k_out/self.k_in)
          self.e_cl = -1*self.RTF*np.log(self.cl_out/self.cl_in)
          #JKKC calculation
          self.j_kcc2 = self.g_kcc2*(self.e_k-self.e_cl)
          #ATPase 
          self.Check_ATP_switch_position()
          self.Calc_ATPase_rate() 
          #Osmolarity and volume adjustments
          self.Calc_new_volume()
          #Incrementing Ion concentration 
          self.Calc_new_ion_conc()
          #Adjustments based based on new volumes
          self.Adjust_concentrations_for_volume()    
          self.volume=self.volume2    
          self.area_scale = self.area/self.volume
          self.f_inv_c_ar = self.F / (self.cm *self.area_scale)

          #Updating Arrays and counters
          if self.t >= self.bin*self.ts :
              self.Update_arrays()
              self.bin += 1       
              
          #Next timestep    
          self.t += self.dt      
