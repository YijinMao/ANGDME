% -------------- The superiority of exploiting sparsity -------------------
%
% This is a simple example to test the Adaptive Nonlinear Group Delay Mode Estimation (ANGDME) algorithm 
%
% Author: Yijin Mao, Hao Liang and Xiaotong Tu
%
% Last modified by: 22/08/20
%

clc; clear; close all

fs = 100;   % sample frequency
T = 15;     % time duration
Nt = 1500;  % the number of samples in time-domain
Nf = floor(Nt/2)+1; % the number of samples in frequency-domain
f = (0 : Nf-1)/T;   % frequency variables
t = (0 : Nt-1)/fs;  % time variables

% Group delay of the signal modes
gd = 0.1*f + 1;

% Amplitude of the signal modes
a = sawtooth(2*pi*f,0.8)+1;

% 'NGDm' denotes a mode of nonlinear group delay signal, 
% 'iNGDFs' denotes its bilateral spectrums, 
% 'ifftSig' denotes its time-domain representation
NGDm = a.*exp(-1j*2*pi*(0.1/2*f.^2 + 1*f + 0.1));  
iNGDFs = [NGDm,conj(fliplr(NGDm(2:ceil(Nt/2))))];  
ifftSig = ifft(iNGDFs);


% time-domain signal
Sig = real(ifftSig);

% frequency-domain signal
NGDmodes = NGDm; 

% initialize
gamma = 1e7;  tol = 2e-4;
iniGD = 3*ones(1,length(f));

% NGDME
tic;
[eGDest, Desest] = ANGDME(NGDmodes,T,iniGD, gamma,tol);
toc;

NGDMEeD = Desest(:,:,end); 


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

%% Reconstructed amplitudes
figure; x1 = 8.5; y1 = -0.2; x2 = 9.5; y2 = 2.2;
plot(f,abs(NGDm),'b','linewidth',2);
hold on;
plot(f,abs(NGDMEeD),'r','linewidth',2);
set(gcf,'Position',[20 100 660 400]);	
set(gcf,'Color','w'); 
xlabel('Frequency/Hz','FontName','Times New Roman');
ylabel('Amplitude','FontName','Times New Roman');
xlim([0,25])
ylim([-1 3]);
set(gca,'FontSize',24)
set(gca,'linewidth',2);
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','k','Linewidth',1);
h1 = axes('position',[0.6 0.3 0.3 0.3]);
axis(h1);
plot(f,abs(NGDm),'b','linewidth',2);
hold on;
plot(f,abs(NGDMEeD),'r','linewidth',2);
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
set(gca,'fontsize',12,'linewidth',1)



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

%% EQF results
Amp1_EQF_NGDME = 20*log10(norm(abs(NGDm) - abs(NGDMEeD),2)/norm(abs(NGDm),2))

GD1_EQF_NGDME = 20*log10(norm(eGDest(:,:,end) - gd,2)/norm(gd,2))

Mode1_EQF_NGDME = 20*log10(norm(real(ifftSig) - real(iffteASig),2)/norm(real(ifftSig),2))


