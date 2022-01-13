%% The main function of this program is to reconstruct the aortic pressure waveform.
close all;
clear all;
clc;
rho=1060;
fss=1000;% New sampling frequency 
fs=128;% Old sampling frequency
% Loading data
Ps=Data_SBP;% Systolic brachial pressure
Pd=Data_DBP;% Diastolic brachial pressure
s=load('Data_pressure.txt');% Aortic and radial pressure waveforms
load('Data_Uwave.mat');% Aortic flow velocity waveform
load('Data_Dmean.mat');% Mean diameter
load('Data_Dmax.mat');% Max diameter
load('Data_Dmin.mat');% Minimum diameter
tau=Data_tau;% Decay constant (can be calculted by the radial pressure waveform)
Dm=dmean/100;% mm to m
Ds=dmax/100;% mm to m
Dd=dmin/100;% mm to m
rP=s(:,3);% Radial pressure waveform
aP=s(:,4);% Aortic pressure waveform
N=length(rP);
t=0:1/fs:1/fs*(N-1);
Pm=mean(rP);
% Cut the descending branch of radial waveform  
nt1=round(0.3*N);
nt2=round(0.45*N);
dn=find(min(rP(nt1:nt2))==rP(nt1:nt2))+nt1-1;
rP_dn=find(rP(dn:end)==max(rP(dn:end)))+nt1-1;
rP_d=rP(rP_dn:end);
t_d=t(rP_dn:end);
Ut=0:1/fss:t(end);
Unt=length(Ut);
U=U(1:Unt);
Umax=max(U);
Utn1=find(U==Umax);
Utn2=find(U(Umax:end)==0);
Utn2=Utn2(2);
% Determining the time of characteristic points for the aortic flow velocitty waveform 
T=t(end);
T1=Ut(Utn1)-Ut(1);
T2=Ut(Utn2)-T1;
T3=Ut(end)-T1-T2;
Ud=U(((T1+T2)*fss+1):end);
t_Ud=(T1+T2):1/fss:T;t_Ud1=0:1/fss:T3;
P0=Pd*exp(T3/tau);
aP_d=P0*exp(-t_Ud1/tau); 
Usn=find(U==max(U));
Us=U(1:Usn);
ts=0:1/fss:(Usn-1)/fss;
% Calculating aortic PWV using Bramwell-Hill method
PPs=(Ds-Dd)*(Pm-Pd)/(Dm-Dd)+Pd;
DD=(Ds^2-Dd^2)/((PPs-Pd)/0.0075*Dd^2);
aPWV=1/sqrt(rho*DD);
% Calculating the waveform at P1 
P1=rho*aPWV*Us*0.0075+Pd; % kpa to mmHg
% Calculating the waveform at P3 
P3=aP_d;
% Calculating the waveform at P2 
M1=mean(P1);
P3=aP_d;
M3=mean(P3);
M2=(Pm*T-M1*T1-M3*T3)/T2;
t_twh=T1;P_twh=P1(end);
t_tn=T1+T2;P_tn=P3(1);
A=[t_twh^2 t_twh 1; t_tn^2 t_tn 1;((t_tn^3)-(t_twh^3))/3 ((t_tn^2)-(t_twh^2))/2 t_tn-t_twh];b=[P_twh; P_tn; M2*T2'];
xx=inv(A)*b; % x=A\b;
a=xx(1);b=xx(2);c=xx(3);
tt=t_twh:1/fss:t_tn;
P2=a*(tt.^2)+b*tt+c;
P_est=[P1(1:end-1) P2(1:end-1) P3];
t_est=0:1/fss:1/fss*(length(P_est)-1);
x=(0:1/fs:1/fs*(N-1))'; 
y = aP;
xi=(0:1/fss:t_est(end))';
yi = interp1q(x,y,xi); 
% Comparison of estimaed and reference pressure waveforms
figure(1)
plot(t_est,P_est)
hold on
plot(t_est,yi,'r')
ylabel('Pressure [mmHg]')
xlabel('Time [s]')
save T T;
% save reAP reAP  % Reference aortic pressure waveform 

