# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:09:13 2020

@author: eshor
"""
class SingleCompartment_PumpLeakModel():
    


    #Class object attribute (i.e. true for all instances of the class)
    R = 8.31466  # J/(Kelvin.mol) #Universal Gas Constant
    F = 96485.0 #C/mol # Faraday Constant in Volts
    T = 310.15   #Kelvin # Absolute temperature (37C)
    RTF = R*T/F  # J/C
    pw = 0.0015 #dm/s #osmotic permeability of biological membranes
    vw = 0.018 #dm^3/mol #partial molar volume of water
    
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    def __init__(self,radius = 5e-5, length=25e-5, capacitance=2e-4): 
    #Attributes
        self.radius = radius
        self.length = length
        self.cm = capacitance
        self.volume = self.np.pi*radius**2*length
        self.area = 2*self.np.pi*radius*length
        self.area_scale = self.area/self.volume
        self.f_inv_c_ar = self.F / (self.cm*self.area_scale)
        
    #Methods
    def Set_InternalConcentrations(self, Na = 14.002e-3, K =122.873e-3 , Cl = 5.163e-3 , impermeant_anions = 154.962e-3,impermeant_anion_charge = -0.85):
        self.na_in = Na
        self.k_in = K
        self.cl_in = Cl
        self.x_in = impermeant_anions
        self.z_x_in = impermeant_anion_charge 
        self.osmol_in = self.na_in+self.k_in+self.cl_in+self.x_in
     
    def Set_ExternalConcentrations(self, Na = 145e-3, K =3.5e-3, Cl = 119e-3,impermeant_anion_charge = -0.85):
        self.na_out = Na
        self.k_out = K
        self.cl_out = Cl
        self.x_out = -1*(self.cl_out - self.na_out - self.k_out)*0.2 
        self.z_x_out = impermeant_anion_charge 
        self.osmol_out = self.na_out + self.k_out +self.cl_out+self.x_out
        
    def Set_Conductances(self, gNa= 2e-3 , gK= 7e-3 ,gCl=2e-3 , gImpermeants=0):
        self.g_na = gNa/self.F
        self.g_k = gK/self.F
        self.g_cl = gCl/self.F
        self.g_x = gImpermeants/self.F
        
    def Set_Timing(self,duration = 10000, time_step = 1e-3, num_intervals = 200):
        self.t = 0                  # real time in seconds
        self.total_t = duration     # total duration
        self.dt = time_step         # time step
        self.total_steps = round(self.total_t/self.dt) #total time steps
        
        self.n = num_intervals      # number of plot points
        self.ctr = 1                # counter for plotting 
        self.ts = self.total_t/self.n  #time in between plot points

    def Set_NaKATPase_Properties(self,time_on=0,time_off=5000,ATP_pump_rate =-1):
        self.atp_t_on = time_on
        self.atp_t_off = time_off
        self.p_default = ATP_pump_rate
        self.p_effective = (10**self.p_default)/self.F
        
    def Set_KCC2_Properties(self,time_on=0,time_off=1000000,KCC2_conductance =2e-3 ):
        self.kcc2_t_on = time_on
        self.kcc2_t_off = time_off
        self.g_kcc2 = KCC2_conductance/self.F
        
    
    def Initialize_Arrays(self):
        self.na_arr = []
        self.k_arr = []
        self.cl_arr = []
        self.x_arr = []
        self.vm_arr = []
        self.t_arr = []           
        self.volume_arr = []
        self.e_cl_arr =[]
        self.e_k_arr = []
            
    def Set_Current(self):
        self.curr = -0*5e-8; #Baseline no current injected

    def Calc_Vm(self):
        self.net_intracellular_charge = (self.na_in+self.k_in)-self.cl_in+self.z_x_in*self.x_in
        self.vm = self.f_inv_c_ar*(self.net_intracellular_charge)
        
    def Simulate(self):
        
        for i in range(1,self.total_steps):  #note tit = total number of timesteps 
        
            # Determining switch position of ATPase
            if (self.t<self.atp_t_off) & (self.t>self.atp_t_on): 
                atp_sw=1 
            else: 
                atp_sw=0
                
            # Voltage Calculation 
            self.Calc_Vm()
            #print(self.vm)
            
            # KCC2 pump rate calculation 
            e_k = self.RTF*self.np.log(self.k_out/self.k_in)
            e_cl = -1*self.RTF*self.np.log(self.cl_out/self.cl_in)
            j_kcc2 = self.g_kcc2*(e_k-e_cl)
            print(e_k)
            print(e_cl)
            print(j_kcc2)
    
            # ATPase pump rate calclation
            j_atp = self.p_effective*(self.na_in/self.na_out)**3
            
            
            #Incrementing Ion concentration    
            d_na = -self.dt*self.area_scale*(self.g_na*(self.vm-self.RTF*self.np.log(self.na_out/self.na_in)) + atp_sw*3*j_atp) 
            d_k = -self.dt*self.area_scale*(self.g_k*(self.vm-self.RTF*self.np.log(self.k_out/self.k_in)) - atp_sw*2*j_atp - j_kcc2) 
            d_cl = self.dt*self.area_scale*(self.g_cl*(self.vm-self.RTF*self.np.log(self.cl_out/self.cl_in)) -j_kcc2) 
            d_x =  self.dt*self.area_scale*(self.g_x*(self.vm-self.RTF*self.np.log(self.x_out/self.x_in)) ) 
            self.na_in += d_na
            self.k_in += d_k
            self.cl_in += d_cl
            self.x_in += d_x
            print(self.na_in)
            print(self.k_in)
            print(self.cl_in)
            print(self.x_in)
            
            #Osmolarity and volume adjustments
            self.osmol_in = self.na_in+self.k_in+self.cl_in+self.x_in
            volume2 = self.volume*(self.osmol_in/self.osmol_out)
    
            #Adjusting Concentrations based on new volumes
            self.na_in *= self.volume/volume2
            self.k_in *= self.volume/volume2
            self.cl_in *= self.volume/volume2
            self.x_in *= self.volume/volume2
            
            self.volume=volume2           
    
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
            print(self.cl_in)
            
    #end of loop
    
    #print statements
        print('Simulation over')
        print('final Cl= ', self.cl_in)

        
        
plm = SingleCompartment_PumpLeakModel()
plm.Set_Conductances()
plm.Set_ExternalConcentrations()
plm.Set_InternalConcentrations()
plm.Set_KCC2_Properties()
plm.Set_NaKATPase_Properties()
plm.Set_Timing()
plm.Set_Current()
plm.Initialize_Arrays()
plm.Calc_Vm()
plm.Simulate()



    
        
"""
########################################################
# 4) SIMULATION #######
########################################################




for i in range(1,totalsteps):  #note tit = total number of timesteps 
    
    
  
            
    # Voltage Calculation 
    Net_Intracellular_Charge = (ConcI_Na+ConcI_K)-ConcI_Cl+z_X*ConcI_X
    Vm = FinvCAr*(Net_Intracellular_Charge)
    
    # KCC2 pump rate calculation 
    EK = RTF*np.log(ConcO_K/ConcI_K)
    ECl = -RTF*np.log(ConcO_Cl/ConcI_Cl)
    JKCC2 = g_KCC2*(EK-ECl)

    # ATPase pump rate calclation
    Jp = p_effective*(ConcI_Na/ConcO_Na)**3

    #Dummy X(impermeant anion) concentrations:
    #ConcO_X_temp = ConcI_X*(1-X_ratio)    

    #Incrementing Ion concentration    
    d_Na = -dt/Ar*(g_Na*(Vm-RTF*np.log(ConcO_Na/ConcI_Na)) + sw*3*Jp) 
    d_K = -dt/Ar*(g_K*(Vm-RTF*np.log(ConcO_K/ConcI_K))-JKCC2-sw*2*Jp) 
    d_Cl = dt/Ar*(g_Cl*(Vm+RTF*np.log(ConcO_Cl/ConcI_Cl))+JKCC2)
    d_X = -dt/Ar*z_X*(g_X*np.log(ConcO_X/ConcI_X))
    ConcI_Na += d_Na   
    ConcI_K += d_K
    ConcI_Cl += d_Cl    
     
    #Osmolarity and volume adjustments
    OsmolI = ConcI_Na+ConcI_K+ConcI_Cl+ConcI_X
    w2 = w*(OsmolI/OsmolO)
    
    #Adjusting Concentrations based on new volumes
    ConcI_Na *= w/w2
    ConcI_K *= w/w2
    ConcI_Cl *= w/w2
    ConcI_X *= w/w2
    
    #Updating Arrays and counters
    w = w2
    
    if t >= ctr*ts :
        Vm_Arr.append(Vm*1e3)
        ECl_Arr.append(ECl*1e3)
        EK_Arr.append(EK*1e3)
        K_Arr.append(ConcI_K*1e3)
        Na_Arr.append(ConcI_Na*1e3)
        Cl_Arr.append(ConcI_Cl*1e3)
        X_Arr.append(ConcI_X*1e3)
        w_Arr.append(100*(1e5)*((3/(4*np.pi))*w)**(1/3))
        t_Arr.append(t)
        ctr += 1
    
        
    t = t+dt      
print('Simulation complete')
print('Final chloride concentration:', Cl_Arr[-1], ' mM')
print('Final cell volume:', w_Arr[-1], 'L')


########################################################

########################################################
# 5) PLOTS
########################################################


plt.plot(t_Arr,Na_Arr,label ="Na")
plt.plot(t_Arr,Cl_Arr,label="Cl")
plt.plot(t_Arr,K_Arr,label="K")
plt.plot(t_Arr,X_Arr,label="X")

SwitchOffX_Arr =[]
SwitchOnX_Arr= []
SwitchOnY_Arr=[]
SwitchOffY_Arr=[]
for a in range(0,150):
  SwitchOffX_Arr.append(t_off)
  SwitchOnX_Arr.append(t_on)
  SwitchOnY_Arr.append(a)
  SwitchOffY_Arr.append(a)

plt.plot(SwitchOnX_Arr,SwitchOnY_Arr,'k:')
plt.plot(SwitchOffX_Arr,SwitchOffY_Arr,'k:')
plt.title("Ionic concentration changes")
plt.xlabel("Time (s)")
plt.ylabel("Intracellular ionic concentrations (M)")
plt.legend(loc = 'upper right')
sns.despine()
plt.annotate('ATPase On', xy=(t_on, 150))
plt.annotate('ATPase Off', xy=(t_off, 150))


plt.plot(t_Arr,EK_Arr,label ="EK")
plt.plot(t_Arr,ECl_Arr,label="ECl")
plt.plot(t_Arr,Vm_Arr,label="Vm")

SwitchOffX_Arr =[]
SwitchOnX_Arr= []
SwitchOnY_Arr=[]
SwitchOffY_Arr=[]
for a in range(-120,150):
  SwitchOffX_Arr.append(t_off)
  SwitchOnX_Arr.append(t_on)
  SwitchOnY_Arr.append(a)
  SwitchOffY_Arr.append(a)

plt.plot(SwitchOnX_Arr,SwitchOnY_Arr,'k:')
plt.plot(SwitchOffX_Arr,SwitchOffY_Arr,'k:')
plt.title("Changes in cell potentials")
plt.xlabel("Time (s)")
plt.ylabel("Voltage (V)")
plt.legend(loc = 'upper right')
sns.despine()
plt.annotate('ATPase On', xy=(t_on, 150))
plt.annotate('ATPase Off', xy=(t_off, 150))


plt.plot(t_Arr,w_Arr,label="Volume (pL)")


SwitchOffX_Arr =[]
SwitchOnX_Arr= []
SwitchOnY_Arr=[]
SwitchOffY_Arr=[]
for a in range(300,530):
  SwitchOffX_Arr.append(t_off)
  SwitchOnX_Arr.append(t_on)
  SwitchOnY_Arr.append(a)
  SwitchOffY_Arr.append(a)
  
plt.plot(SwitchOnX_Arr,SwitchOnY_Arr,'k:')
plt.plot(SwitchOffX_Arr,SwitchOffY_Arr,'k:')
sns.despine()

plt.annotate('ATPase On', xy=(t_on, max(w_Arr)))
plt.annotate('ATPase Off', xy=(t_off, max(w_Arr)))"""
