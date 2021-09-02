%% PLM_CD- Pump-Leak model - CD integration of Keener-Sneyd model
% p=pump rate., X = moles of impermeant ion, z = its charge 
% R=RT/F, all concens in M, V in volts, dimensions dm. 
% start at area & volume defined by radius, the volume is allowed to change
% but the area is fixed 
%%%%%%%%%%%%%%%%%%%

clear, clc; 
R=25.69*1e-3; F=96485; % R (RT/F) in Volts, Faraday's constant C/mol
n=200; % # points to plot 
Vm=1:n-1; K=1:n-1; Na=1:n-1; Cl=1:n-1; W=1:n-1; X=1:n-1; time=1:n-1;%create plotting arrays
time(1)=0; dt=1e-3;  %zero time, dt time step, n number of points to plot
tt=10000; ts=tt/n; tit=round(tt/dt);% tt total time, tit - total # of steps 
ton=0; toff=5000; %time when pump turned on & off 
curr=-0*5e-8; % current injected
gna=0.01*0.1/F ; gk=0.3*0.1/F ; gcl=0.2*0.1/F; % conductances in mS/cm^2 conv to S/dm^2 (10^-3/10^-2) 
p=0.5e-4/F; %pump rate (C/dm^2 s) div by F 
ck=2;cna=3;%pump stoichiometry
rad=5; rad=rad*1e-5; %(rad in um convert to dm)
w=(4/3)*pi*rad^3; %vol in liters
Ar=4*pi*rad^2; %area in dm^2
C=1e-4; %capacitance F/dm^2
FinvCAr=F/(C*Ar); %(F/C*Area)
Oso=300e-3; ko=3e-3; % set Oso 
xo=0, zo=0; %xo is fixed ion outside and zo its charge
clo=(Oso-xo+zo*xo)/2;% calc clo from charge and osmotic balance 
nao=clo-ko-zo*xo; % calc nao 
x=50e-3; z=-1; MX=x*w;% X CONCEN
cl=(Oso-2*x)/2 ,  na=0.8*(Oso-cl-x), k=0.2*(Oso-cl-x), %INITIAL INTRA CONCs
ctr=1;t=0; % ctr counter for plotting points, t real time 5%%
sw=0;
x=MX/w;
V =FinvCAr*w*(na+k-cl+z*x) % starting voltage
for i = 2:tit 
    if (toff>t)&&(t>ton) sw=1 ;
else sw=0;
   if t>toff sw=0;      
    end; 
end;
V =FinvCAr*w*(na+k-cl+z*x);
invw=1/w;
dna=-dt*Ar*invw*(gna*(V-R*log(nao/na))+sw*cna*p);%increase in Na during time step dt
dk=-dt*Ar*invw*(gk*(V-R*log(ko/k))-sw*ck*p+sw*curr);
dcl=dt*Ar*invw*(gcl*(V+R*log(clo/cl)));
na=na+dna; k=k+dk; cl=cl+dcl; % increment concentrations

  Osi=na+k+cl+x; % intracellular osmolarity 
 w2=(w*Osi)/Oso; % correct volume 
na=(na*w)/w2; k=(k*w)/w2; cl=(cl*w)/w2 ; x=(x*w)/w2; % correct concens 
w=w2;

  if t>=ctr*ts 
Vm(ctr)=1000*V; K(ctr)=1000*k; Na(ctr)=1000*na; 
Cl(ctr)=1000*cl; W(ctr)=100*(1*10^5)*((3/(4*pi)).*w).^(1/3);X(ctr)=1000*x; time(ctr)=t; ctr=ctr+1;
 end
 t=t+dt;
 end
figure
subplot(2,1,1)
plot (time,K,'--k',time,Na,'-r', time,Cl,'-g',time,X,'-b')
legend('K','Na','Cl','X')
subplot(2,1,2);
plot(time,W);
legend('Rad(um x100)','Vm')
