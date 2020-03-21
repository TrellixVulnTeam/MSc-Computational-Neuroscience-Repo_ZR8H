# -*- coding: utf-8 -*-
"""
Created on Sun Aug 14 15:43:20 2016

@author: Kira

PLM model including KCC2

:::FUNCTIONS:::

plm(p,graph,...)
==> runs a time series plm run
p := desired pump rate
graph := {1 iff a graph is desired as output; and any other value otherwise}
various other parameters to play with - backbone of standard simulations
time (dt / t) is specified in seconds

zplm(z,...) 
==> runs the parametric solution over log pump rates in P for impermeant anion charge of z
z := charge of impermeant anions
gk,gna,gcl,gkcc := (conductances) can set desired conductance values
molinit := initial total mols in the cell (determines volume)

zp(Z,p,...) 
==> runs the parametric solution at the pump rate defined by p for different z values
Z := desired array of impermeant anion charges, multiplied by 100 (i.e. input as range(start*100, end*100))
p := initial log pump value satisfying the pump rate given by 10**(p)/(F*R)
gk,gna,gcl,gkcc := (conductances) can set desired conductance values
molinit := initial total mols in the cell (determines volume)
moldelt := set how much total mols change by across the simulation

ecl(cli)
==> returns Nernst chloride potential
"""

import numpy as np

from pylab import rcParams
rcParams['figure.figsize'] = 8,8


# constants, fixed parameters
R=26.725*1e-3 #E: R in this context is not the gas constant, it's RT/F
F=96485.0 # R (RT/F) in Volts, where F is Faraday's constant in C/mol, and T is 37 deg C
n=200 # points to plot 
gna=2e-3/F #E: why are conductances/F?
gk=7e-3/F
gcl=2e-3/F # gna,gk,gcl: conductances in mS/cm^2 conv to S/dm^2 (10^-3/10^-2) - corrected for neuron
gkcc=2e-3/F # gkcc conductance
ck=2
cna=3 # cna,ck: pump (ATPase) stoichiometries
rad=5*1e-5 # radius in um convert to dm
rad0=rad
length=25*1e-5 # length in um converted to dm
nao=145e-3
clo=119e-3
ko=3.5e-3 # nao,clo,ko: extracellular concentrations (mM converted to M)
z=-0.85 # intracellular (and extracellular) charge of impermeant anions
gamma=gna/gk
beta=1.0/(gk*gcl+gkcc*gk+gcl*gkcc)
nae=nao
ke=ko
cle=clo
xe1=-1*(cle-nae-ke)
xe=xe1*0.2
ose=xe1+cle+nae+ke # extracellular osmolarity
P=range(-70000,-38000)
default_p=-1
default_P=-40456 # P_effective x10^5
vw=0.018 # partial molar volume of water, dm3/mol
pw=0.0015 # osmotic permeability, biological membrane, dm s
km=6*10**(-7) # extensional rigidity of RBC at 23 deg, N/dm
km2=2.5*10**(1)
density=1.0 # kg/dm3 = g/ml --> assume close to 1 (density of water)
hp=1e-3
hydrop=0
qpump=6.13*1e-5 #picoamperes
kd=15*1e-3 #M Kd (Raimondo 2012)
vmax=5*1e-3 #M/s Vmax (Raimondo 2012)

w=np.pi*rad**2*length # initial volume in liters
sa=2*np.pi*rad*(length)
w1=w # initial volume stored for graphing later
Ar=2.0/rad # area constant (F and H method)
C=2e-4 # capacitance (F/dm^2)
FinvCAr=F/(C*Ar) # (F/C*area scaling constant)
V=FinvCAr*(14.002e-3+122.873e-3 - 5.16e-3 - 154.962e-3)
"""
def plm(p=(10**(default_p))/(F),graph=0,pkcc=gkcc,gx=0,xt=100000,os_init=ose,clinit=5.163e-3,toff=150000,ton=150000,tt=200,xinit=154.962e-3,two=0,xe=xe,f4d=0,ke=ke,n=1800,k_init=122.873e-3,na_init=14.002e-3,tk=100000,ratio=0.98,xend=120,osmofix=False,paratwo=False,moldelt=1e-13,xflux=0,z=z,dz=0,Zx=-1,ztarget=-100,length=length,areascale=1,rad=rad,title='fig.eps',neww=0,ls='-',a0=0,a1=0,a2=0,os_choose=0,f1d=False,hamada=0,kccmodel=0,vmax=vmax,lin=0):
    # create plotting arrays
    Vm=[]
    K=[]
    Na=[]
    Cl=[]
    W=[]
    X=[]
    time=[]
    Cl2=[]
    X2=[]
    Na2=[]
    K2=[]
    z_delt=[]
    xe_delt=[]
    gkcc_delt=[]
    kflux=[]
    naflux=[]
    clflux=[]
    Xflux=[]
    wflux=[]
    
    dt=1e-3 # zero time, dt time step
    ts=tt/n # plotting timestep 
    ctr=1 # counter for plotting points
    t=0 # real time
    sw=1 # switch for ATPase action 
    
    
    if f1d==True:
        w=w*154.962e-3/xinit # adjust for starting conditions in F1D (optimisation)
    
    if areascale==0 or areascale==1:
        Ar=sa/w
   
    sarest=sa
    
    # na,k,cl,x: intracellular starting concentrations
    na=na_init
    x=xinit
    #cl=((os_init-na-k)*z+na+k)/(1+z)
    cl=clinit
    #x=(cl-na-k)/z
    k=k_init
    cle=clo
    
    dw=0
    dk=0
    dcl=0
    dna=0
    dx=0
        
    if osmofix==True:
        if xinit==0:
            x=(os_init-2*cl)/(1-z)
            xinit=x
        else:
            cl=(os_init+(z-1)*x)/2.0
            print(cl)
    
    if k_init==0:
        k=cl-z*x-na
        
    print("k_init: "+str(k))
    print("osi: "+str(k+cl+x+na))
    print("z_aim: "+str(ztarget) +" with zflux of "+str(Zx))
    
    xm=x*ratio
    xtemp=x*(1-ratio)
    zxm=z
    zx=z
    
    # for f1c --> slow change in ATPase rate
    pdinit=-5.0
    pd=default_p
    em=(default_p-pdinit)/(12.0*10**4)/8
    jeffconstant=p*(na/nao)**3
    
    # related to anion flux
    if two==1:
        zx=Zx
        zxm=(z*x-zx*xtemp)/xm
        if paratwo==True:
            return (w*x)
    
    while t < tt: # loop over time              
        V=FinvCAr*(na+k-cl+z*x) # voltage
        
        # update arrays for plotting
        if t>=(ctr-1)*ts:
            K.append(1000*R*np.log(ke/k))
            K2.append(1000*k)
            Na.append(1000*R*np.log(nao/na))
            Na2.append(1000*na)
            Cl.append(1000*R*np.log(cl/cle))
            Cl2.append(cl*1000.0)
            X.append(z*1000*R*np.log(xe/x))
            X2.append(1000*(x))
            W.append(w)
            Vm.append(1000*V)
            time.append(t)
            z_delt.append(z)
            xe_delt.append(xe)
            gkcc_delt.append(pkcc)
            wflux.append(dw)
            Xflux.append(dx)
            clflux.append(dcl)
            naflux.append(dna)
            kflux.append(dk)
            ctr+=1
        
        # various conditional states
        if tk+360>t>tk:
            pkcc+=1e-12    # control switch for gkkc ramp (Fig 3)
            vmax+=3.3e-7
        
        if dz!=0 and xt<t<xt+420 and xtemp>0 and xm>0:
            xtemp+=dz
            xm-=dz # control switch for anion flux by modifying z directly

        if two==1:
            z=(zxm*xm+zx*xtemp)/(xm+xtemp) # recalculate average charge if needed
        
        if f4d!=0:
            if xt+400>t>xt:
                xe+=f4d*6e-5
                cle-=f4d*6e-5 # Figure 4D (balance the charge differences) --> can adjust the ratio at * for interest
        
        jp=p*(na/nao)**3 # cubic pump rate update (dependent on sodium gradient)
        
        if hamada!=0:
            jp=qpump*hamada*(1.62/(1+(0.0067/na)**3)+1.0/(1+(0.0676/na)**3))/F
        
        if lin!=0:
            jp=p*(na/nao)
        
        if neww==4 or neww==5:
            jp=jeffconstant # Figure 6
            
        if (toff>t) and (t>ton):
            if pd>pdinit:
                pd-=em
                p=(10**(pd))/F
        elif t>toff:
            if pd<default_p:
                pd+=em
                p=(10**(pd))/F # ATPase ramp

        # kcc2
        if kccmodel==1: 
            externals=cle*ke/k
            fix=np.log(externals/0.056)
            jkcc2=-0.117*np.log(externals/cl)/fix*vmax*cl/(kd+cl)/F #Raimondo
        elif kccmodel==2:
            jkcc2=51.55*pkcc*(ke*cle-k*cl) #Fraser and Huang
        elif kccmodel==3:
            jkcc2=0.011125*0.3*(ke*cle-k*cl)/(0.000054*((1+ke*cle/0.000054)*(1+ke/0.009)*(1+cle/0.006)+(1+k*cl/0.000054)*(1+k/0.009)*(1+cl/0.006)))/F
        else:
            jkcc2=pkcc*(K[ctr-2]-Cl[ctr-2])/1000.0 #Doyon

        # ionic flux equations
        dna=-dt*Ar*(gna*(V-R*np.log(nao/na))+cna*jp*sw) 
        dk=-dt*Ar*(gk*(V-R*np.log(ke/k))-ck*jp*sw-jkcc2)
        dcl=dt*Ar*(gcl*(V+R*np.log(cle/cl))+jkcc2) #dna,dk,dcl: increase in intracellular ion conc during time step dt
        dx=-dt*Ar*zx*(gx*(V-R/zx*np.log(xe/(xtemp))))
        na+=dna
        k+=dk
        cl+=dcl # increment concentrations
        
        # anion flux switches
        if xend==0 and (t>xt):
            if (np.abs(x*w-xinit*w1)<moldelt) and (abs((np.abs(z)-np.abs(ztarget)))>0.001) and (min(z,zx)<=ztarget<=max(z,zx)):
                if xflux==0:
                    xtemp+=dx
                    tt=t+180
                else:
                    xtemp+=xflux
                    dx=xflux
                    tt=t+1000
            else:
                if (min(z,zx)<=ztarget<=max(z,zx)):
                    print('anions stopped diffusing at '+str(t))
                    xend=1
                    dx=0
                else:
                    if xflux!=0 and tt<1000 and xtemp>0 and (min(zxm,zx)<=ztarget<=max(zxm,zx)):
                        xtemp-=xflux
                        tt=t+50
                        dx=-xflux
                    else:
                        print( 'anions stopped diffusing at '+str(t))
                        xend=1
                        dx=0
                
        if xt+xend>t>xt:
            if xflux!=0:
                xtemp+=xflux/10
                dx=xflux/10
            else:
                xtemp+=dx/10 
                dx=dx/10
        elif xend!=0:
            dx=0
        elif t<xt:
            dx=0
            
        # update volume (usual method)
        x=xm+xtemp
        osi=na+k+cl+x # intracellular osmolarity 
        ose=nae+ke+cle+xe+xe1*0.8
        dw=dt*(vw*pw*sa*(osi-ose))
        w2=w+dw
        
        # other volume updates (incorporating hydrostatic pressure - various options considered) 
        if neww==1:
            w2=w+dt*(vw*pw*sa*(osi-ose)+hp*dt/density*km*(sarest-sa)/sarest)
        elif neww==2:
            w2=w+dt*(vw*pw*sa*(osi-ose-os_choose))
        elif neww==3 or neww==5:
            hydrop=4.0*km2*np.pi*(rad/rad0-1)/(R*F)
            w2=w+dt*(vw*pw*sa*(osi-ose-hydrop))

        # correct ionic concentrations and surface area by volume change
        na=(na*w)/w2
        k=(k*w)/w2
        cl=(cl*w)/w2
        x=(x*w)/w2
        xm=(xm*w)/w2
        xtemp=(xtemp*w)/w2
        w=w2
        sa=2*np.pi*rad*(length)
        
        # methods of updating Ar constant (dependent on how the surface area changes for volume, by radius or length)
        if areascale==1:
            rad=np.sqrt(w/(np.pi*length))
            Ar=sa/w
            FinvCAr=F/(C*Ar)
        elif areascale==0:
            length=w/(np.pi*rad**2)
            
        t+=dt

    print('na', na, 'k', k, 'cl', cl, 'x', x, 'vm', V, 'cle', cle, 'ose', ose, 'osi', osi, 'deltx', x*w-xinit*w1)
    print('w', w, 'radius', rad, 'z', z)
    print('ecl', Cl[-1])
    return na, k, cl, x, V, Na[-1], K[-1], Cl[-1], X[-1], Vm[-1], W, time, Na, K, Cl, X, Vm, Cl2, Na2, K2, X2, w, z_delt, xe_delt, gkcc_delt, a0, a1, a2, naflux, kflux, clflux, wflux, Xflux, np.log10(jp*F), osi, ose

"""
