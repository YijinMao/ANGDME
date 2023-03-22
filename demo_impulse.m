% ----------------- Transient impulse signal analysis ---------------------
%
% This is a simple example to test the Nonlinear Group Delay Mode Estimation (NGDME) algorithm 
% 
% Author: Hao Liang and Xiaotong Tu
%
% Last modified by: 21/08/23
%

clc; clear; close all

fs = 512;   % sample frequency
Nt = 512 + 1; % the number of samples in time-domain
t = (0:Nt-1)/fs;  %time variables
f = (0:(Nt-1)/2)*fs/Nt; % frequency variables
w = 2*pi*f; % auxiliary variable
T = t(end); % time duration

% Transient impulse signal 1
T_1 = 0.4; T_2 = (0.4-0.2)/(f(end)-f(1));
T_3 = 5; T_4 = -0.001; T_5 = 0.01;
NGDm = (20*exp(-0.5/800000*(w-100*2*pi).^2)).*exp(-1j*(T_1*w+1/2*(T_2/(2*pi))*w.^2+T_3*exp(T_4*w).*sin(T_5*w)));
gd = T_1+T_2*f+T_3*exp(T_4*2*pi*f).*(T_5*cos(T_5*2*pi*f)+T_4*sin(T_5*2*pi*f));
iNGDFs = [NGDm conj(fliplr(NGDm(2:end)))]; Sig = ifft(iNGDFs);


%% NGDME
gamma = 1e4; lambda = 2e-2; tol = 1e-8;
iniGD = 0.48*ones(1,length(f));

tic;
[eGDest, Desest] = ANGDME(NGDm,T,iniGD,gamma,tol);
toc;

NGDMEeD = Desest(1,:,end); 

% Obtain the time-domain signal by inverse FFT
ieADFs = [NGDMEeD,conj(fliplr(NGDMEeD(2:ceil(Nt/2))))]; iffteASig = ifft(ieADFs); 

%% initial group delay
figure
plot(iniGD,f,'b','linewidth',3) % initial group delay
set(gcf,'Position',[20 100 640 500]);	 
xlabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
ylabel('Time/s','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'linewidth',2);
set(gcf,'Color','w');
axis([0 1 0 250])
%% Reconstructed modes
figure
set(gcf,'Position',[20 100 640 500]);	     
set(gcf,'Color','w'); 
plot(t,real(Sig),'b','linewidth',2);   % time-domain signal modes
hold on
plot(t,real(iffteASig),'r--','linewidth',2);
xlabel('Time/s','FontSize',24,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',24,'FontName','Times New Roman');
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);
ylim([-6,6])

%% Estimated group delay by NGDME
figure
plot(gd,f,'b','linewidth',3) % true group delay
hold on
plot(eGDest(:,:,end),f,'r--','linewidth',3) % estimated group delay
set(gcf,'Position',[20 100 640 500]);	 
xlabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
ylabel('Time/s','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'linewidth',2);
set(gcf,'Color','w');	
axis([0 1 0 250])

%% EQF results
Mode_EQF_NGDME = 20*log10(norm(real(Sig) - real(iffteASig),2)/norm(real(Sig),2))
GD_EQF_GDMD = 20*log10(norm(eGDest(:,:,end) - gd,2)/norm(gd,2))


