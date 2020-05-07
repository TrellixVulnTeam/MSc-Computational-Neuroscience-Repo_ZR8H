# -*- coding: utf-8 -*-
"""
Created on Wed Apr 22 16:05:53 2020

@author: eshor
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class SingleCompPLM():

  R = 8.31466  # J/(Kelvin.mol) #Universal Gas Constant
  F = 96485.0 #C/mol # Faraday Constant in Volts
  T = 310.15   #Kelvin # Absolute temperature (37C)
  RTF = R*T/F  # J/C
  pw = 0.018 #dm/s #Osmotic permeability
  vw = 0.018 #dm3/mol #Partial molar volume of water
  default_p = 0.1 #C/dm2.s #default pump rate

  def __init__(self, radius = 5e-5, length= 25e-5, capacitance = 2e-4):
      self.radius = radius
      self.length = length
      self.cm = capacitance
      self.volume = np.pi*self.radius**2*self.length
      self.area = 2*np.pi*self.radius*self.length
      self.area_scale = self.area/self.volume
      self.f_inv_c_ar = self.F / (self.cm *self.area_scale)


  def SetConcentrations(self, Na = [14.002e-3, 145e-3], K = [122.873e-3, 3.5e-3], Cl = [5.163e-3, 119e-3], Impermeants = [154.962e-3, 29.5e-3], Impermeant_charges = [-0.85,-0.85] ):
      self.na_in = Na[0]
      self.na_out = Na[1]
      self.k_in = K[0]
      self.k_out = K[1]
      self.cl_in = Cl[0]
      self.cl_out = Cl[1]
      self.x_in = Impermeants[0]
      #self.x_out = Impermeants[1]
      self.x_out = -1*0.2*(self.cl_out-self.na_out-self.k_out)
      self.z_in = Impermeant_charges[0]
      self.z_out = Impermeant_charges[1]
      # Kira uses these to get starting extravellular impearmeants xe1=-1*(cle-nae-ke)
      # xe=xe1*0.2


  def SetConductances(self, gNa= 2e-3 , gK= 7e-3 ,gCl=2e-3 , gImpermeants=0):
      self.g_na = gNa/self.F
      self.g_k = gK/self.F
      self.g_cl = gCl/self.F
      self.g_x = gImpermeants/self.F

  def Initialize_Arrays(self):
      self.na_arr = []
      self.k_arr = []
      self.cl_arr = []
      self.x_arr = []
      self.vm_arr = []
      self.vm_test = []
      self.na_test =[]
      self.k_test = []
      self.cl_test = []
      self.x_test =[]
      self.t_arr = []           
      self.volume_arr = []
      self.e_cl_arr =[]
      self.e_k_arr = []


  def CalcVm(self):
       na = self.na_in
       k = self.k_in
       cl = self.cl_in
       z = self.z_in
       x = self.x_in
       xz = x*z
       net_charge = float(na + k - cl + (xz))
       self.vm = self.f_inv_c_ar*(net_charge)

  def Set_NaKATPase_Properties(self,time_off=300,time_on=900,ATP_pump_rate =-1):
       self.atp_t_on = time_on
       self.atp_t_off = time_off
       self.p_default = ATP_pump_rate
       self.p_effective = (10**self.p_default)/self.F
       
 
  def Set_KCC2_Properties(self,time_on=0,time_off=1000000,KCC2_conductance =2e-3 ):
       self.kcc2_t_on = time_on
       self.kcc2_t_off = time_off
       self.g_kcc2 = KCC2_conductance/self.F

  def SetTiming(self,duration = 1000, time_step = 1e-3, num_intervals = 180000):
       self.t = 0                  # real time in seconds
       self.total_t = duration     # total duration
       self.dt = time_step         # time step
       self.total_steps = round(self.total_t/self.dt) #total time steps
       self.n = num_intervals      # number of plot points
       self.ctr = 1                # counter for plotting 
       self.ts = self.total_t/self.n  #time in between plot points   
       
   
       
  def Simulate(self):
      
      d_na =0
      d_k =0
      d_cl =0
      volume2=0
    
      for i in range(1, self.total_steps):  
          
          if (self.t>self.atp_t_off) & (self.t<self.atp_t_on):
                atp_sw = 0
          else:
                atp_sw = 1     
            
          self.CalcVm()
          self.vm_test.append(self.vm)
          e_k = self.RTF*np.log(self.k_out/self.k_in)
          e_cl = -1*self.RTF*np.log(self.cl_out/self.cl_in)
          
          self.na_test.append(self.na_in)
          self.k_test.append(self.k_in)
          self.cl_test.append(self.cl_in)
          self.x_test.append(self.x_in)
          
          #JKKC calculation
          j_kcc2 = self.g_kcc2*(e_k-e_cl)
          
          
          #ATPase ramp
          pdinit = -5 
          pd = self.p_default
          ramp = (pd-pdinit)/(12*10**4)/8   #equivalent to em in Kira's code
          
          if atp_sw  ==1:
              if pd>pdinit:
                  pd-=ramp
                  self.p_effective = (10**pd)/self.F
          else:
              if pd<pdinit:
                  pd += ramp
                  self.p_effective = (10**pd)/self.F
                  
          j_atp = self.p_effective*(self.na_in/self.na_out)**3
          
          #Osmolarity and volume adjustments
          self.osmol_in = self.na_in+self.k_in+self.cl_in+self.x_in  
          self.osmol_out = self.na_out + self.k_out +self.cl_out+self.x_out 
          self.radius = np.sqrt(self.volume/(np.pi*self.length))
          self.area = 2*np.pi*self.radius*self.length 
          dw=self.dt*self.vw*self.pw*self.area*(self.osmol_in-self.osmol_out)  
          volume2 = self.volume+dw
          
          
          #Incrementing Ion concentration    
          d_na = -self.dt*self.area_scale*(self.g_na*(self.vm-self.RTF*np.log(self.na_out/self.na_in)) + atp_sw*3*j_atp) # - (1/self.volume)*dw*self.na_in
          d_k = -self.dt*self.area_scale*(self.g_k*(self.vm-self.RTF*np.log(self.k_out/self.k_in)) - atp_sw*2*j_atp - j_kcc2)#  - (1/self.volume)*dw*self.k_in
          d_cl = self.dt*self.area_scale*(self.g_cl*(self.vm-self.RTF*np.log(self.cl_in/self.cl_out)) - j_kcc2) # - (1/self.volume)*dw*self.cl_in
          self.na_in += d_na
          self.k_in += d_k
          self.cl_in += d_cl
                 
          
          #Adjusting Concentrations based on new volumes
          self.na_in *= self.volume/volume2
          self.k_in *= self.volume/volume2
          self.cl_in *= self.volume/volume2
          self.x_in *= self.volume/volume2
                  
          self.volume=volume2    
           
          
          self.area_scale = self.area/self.volume
          self.f_inv_c_ar = self.F / (self.cm *self.area_scale)

        
          #Updating Arrays and counters
          if self.t >= self.ctr*self.ts :
              self.vm_arr.append(self.vm*1e3)
              self.e_cl_arr.append(e_cl*1e3)
              self.e_k_arr.append(e_k*1e3)
              self.na_arr.append(self.na_in*1e3)
              self.k_arr.append(self.k_in*1e3)
              self.cl_arr.append(self.cl_in*1e3)
              self.x_arr.append(self.x_in*1e3)
              self.volume_arr.append(self.volume)
              self.t_arr.append(self.t)
              self.ctr += 1       
              
                    
          self.t += self.dt      

#####################
plm = SingleCompPLM()
plm.SetConcentrations()
plm.SetConductances()
plm.SetTiming()
plm.Set_KCC2_Properties()
plm.Set_NaKATPase_Properties()
plm.Initialize_Arrays()
plm.Simulate()

"""
plt.plot(plm.t_arr,plm.na_arr,label ="Na")
plt.plot(plm.t_arr,plm.cl_arr,label="Cl")
plt.plot(plm.t_arr,plm.k_arr,label="K")
plt.plot(plm.t_arr,plm.x_arr,label="X")
"""

plt.plot(plm.t_arr,plm.vm_arr,label ="Na")
"""
plt.plot(plm.t_arr,plm.e_cl_arr,label="Cl")
plt.plot(plm.t_arr,plm.e_k_arr,label="K")
"""