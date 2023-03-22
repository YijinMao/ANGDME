% --------------------- Cross-mode signal analysis ------------------------
%
% This is a simple example to test the Nonlinear Group Delay Mode Estimation (NGDME) algorithm 
% 
% Author: Hao Liang and Xiaotong Tu
%
% Last modified by: 21/08/23
%

clc; clear; close all

fs = 100;   % sample frequency
T = 15;     % time duration
Nt = 1500;  % the number of samples in time-domain
Nf = floor(Nt/2)+1; % the number of samples in frequency-domain
f = (0 : Nf-1)/T;   % frequency variables
t = (0 : Nt-1)/fs;  % time variables

% Group delay of the signal modes
gd1 = -0.001*f.^2 + 0.3*f+1;
gd2 = -1/250*f.^2  + 12;

% Amplitude of the signal modes
a1 = (1 + 0.2*sin(2*pi*2/50*f)); 
a2 = (1 + 0.2*cos(2*pi*2/50*f)); 

% 'NGDm1' denotes a mode of nonlinear group delay signal, 'iNGDFs1' denotes 
% its bilateral spectrums, and 'ifftSig1' denotes its time-domain representation
NGDm1 = a1.*exp(-1j*2*pi*(-0.001/3*f.^3 + 0.3/2*f.^2 + 1*f + 0.1));  
iNGDFs1 = [NGDm1,conj(fliplr(NGDm1(2:ceil(Nt/2))))]; ifftSig1 = ifft(iNGDFs1);
NGDm2 = a2.*exp(-1j*2*pi*(-1/750*f.^3 + 12*f + 0.8)); 
iNGDFs2 = [NGDm2,conj(fliplr(NGDm2(2:ceil(Nt/2))))]; ifftSig2 = ifft(iNGDFs2);


% time-domain signal
Sig1 = real(ifftSig1);
Sig2 = real(ifftSig2);
Sig = Sig1+Sig2;

% frequency-domain signal
NGDmodes = NGDm1 + NGDm2; 

%% TSST1
stft_type ='T';
tf_type ='STFT';
direction = 'T';
gamma = 0.00001;
tradeoff = 0.009;
tf_type ='TSST1';s = 5e-1;
[~,TFx,~,~,Rep,~,q,~,~] = TFM(Sig,fs,s,stft_type,tf_type,gamma);
TSET1_TFD = abs(TFx');
figure; x1=7.5; y1=23; x2=10; y2=33;
set(gcf,'Position',[20 100 640 500]);
imagesc(t,f,TSET1_TFD);
axis([0 15 0 50]);
xlabel('Time/s','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'YDir','normal')
%%zoom area
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',1);
h1=axes('position',[0.6 0.2 0.3 0.3]);
axis(h1);
imagesc(t,f,TSET1_TFD);
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]); %%不显示刻度
set(gca,'fontsize',12,'linewidth',1,'XColor','r','YColor','r') %%XY坐标轴的颜色
colormap turbo %% 颜色图的色阶
%% STFT
stft_type ='T';
tf_type ='STFT';
direction = 'T';
gamma = 0.00001;
tradeoff = 0.009;

s = 2e-1;
[~,TFx,~,~,Rep,~,q,~,~] = TFM(Sig,fs,s,stft_type,tf_type,gamma);
STFT_TFD = abs(TFx');
figure; x1=7.5; y1=23; x2=10; y2=33;
set(gcf,'Position',[20 100 640 500]);
imagesc(t,f,STFT_TFD);
axis xy
xlabel('Time/s','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency/Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'YDir','normal')
%%zoom area
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',1);
h1=axes('position',[0.6 0.2 0.3 0.3]);
axis(h1);
imagesc(t,f,STFT_TFD);
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]); %%不显示刻度
set(gca,'fontsize',12,'linewidth',1,'XColor','r','YColor','r') %%XY坐标轴的颜色
colormap turbo %% 颜色图的色阶
%% ridge extraction
bw = T/100;   % the bandwidth of the TF filter for ridge extraction
beta1 = 1e-4; % beta1 should be larger than the following beta
num = 2;      % the number of the components
delta = 20; Nfrebin = 1024; window = 64;
[tidexmult, tfdv] = extridge_mult(NGDmodes, fs, num, delta, beta1, bw, Nfrebin,window);

%% ridge path regrouping (RPRG)
% the RPRG algorithm is developed to extract ridge curves of crossed signal modes
% more details about the RPRG can be found in paper: Chen S, Dong X, Xing G, et al, Separation of Overlapped Non-Stationary Signals by Ridge Path Regrouping and Intrinsic Chirp Component Decomposition, IEEE Sensors Journal, 2017.
thrt = length(t)/30;
[tindex,~] = RPRG(tidexmult,thrt);

%% NGDME
gamma = 1e7;  tol = 1e-8;

iniGD = curvesmooth(t(tindex),1e-7); % smooth the GD curve
% iniGD = [iniGD(1,:)+2; iniGD(2,:)];
tic;
[eGDest, Desest] = ANGDME(NGDmodes,T,iniGD,gamma,tol);
toc;

NGDMEeD1 = Desest(1,:,end); NGDMEeD2 = Desest(2,:,end); 

% Obtain the time-domain signal by inverse FFT
ieADFs1 = [NGDMEeD1,conj(fliplr(NGDMEeD1(2:ceil(Nt/2))))]; iffteASig1 = ifft(ieADFs1); 
ieADFs2 = [NGDMEeD2,conj(fliplr(NGDMEeD2(2:ceil(Nt/2))))]; iffteASig2 = ifft(ieADFs2);
eSig1 = real(iffteASig1);
eSig2 = real(iffteASig2);

%% EQF results
Mode1_EQF_NGDME = 20*log10(norm(real(ifftSig1) - real(iffteASig1),2)/norm(real(ifftSig1),2))
Mode2_EQF_NGDME = 20*log10(norm(real(ifftSig2) - real(iffteASig2),2)/norm(real(ifftSig2),2))


%% FIGURE %%
%% time-domain signal
figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,Sig,'b','linewidth',2);
xlabel('Time (s)','FontSize',24,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
set(gca,'FontName','Times New Roman')
axis([0 15 -0.2 0.2]);

%% initial group delay
figure
plot([gd1; gd2],f,'b','linewidth',3) % initial group delay
set(gcf,'Position',[20 100 640 500]);	 
ylabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
xlabel('Time/s','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'linewidth',2);
set(gcf,'Color','w');
axis([0 15 0 50])

%% Reconstructed group delays
figure; 
% x1=7.5; y1=23; x2=10; y2=33;
x1=11; y1=0.5; x2=13.5; y2=12;
plot([gd1;gd2],f,'b','linewidth',2);
hold on;
plot(eGDest(:,:,end),f,'r','linewidth',2);
set(gcf,'Position',[20 100 640 500]);	
set(gcf,'Color','w'); 
xlabel('Time/s','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
ylim([0 50]);
legend('Original','', 'Estimated')
set(gca,'FontSize',24);
set(gca,'linewidth',2);
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','k','Linewidth',1);
h1=axes('position',[0.6 0.2 0.3 0.3]);
axis(h1);
plot([gd1;gd2],f,'b','linewidth',1.5);
hold on;
plot(eGDest(:,:,end),f,'r','linewidth',1.5);
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]);
set(gca,'fontsize',12,'linewidth',1)

%% Time-frequency representation
%% ANGDME 
Ampest = abs(Desest(:,:,end));
timerange = [0 T];
[ASpec tbin] = TFspec(eGDest,Ampest,timerange);
figure; x1=7.5; y1=23; x2=10; y2=33;
set(gcf,'Position',[20 100 640 500]);	
imagesc(tbin,f,abs(ASpec)); 
axis xy
xlabel('Time/s','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency/Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'YDir','normal')
%%zoom area
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',1);
h1=axes('position',[0.6 0.2 0.3 0.3]);
axis(h1);
imagesc(t,f,abs(ASpec));
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]); %%不显示刻度
set(gca,'fontsize',12,'linewidth',1,'XColor','r','YColor','r') %%XY坐标轴的颜色
colormap turbo %% 颜色图的色阶

%% STFT
stft_type ='T';
tf_type ='STFT';
direction = 'T';
gamma = 0.00001;
tradeoff = 0.009;

s = 2e-1;
[~,TFx,~,~,Rep,~,q,~,~] = TFM(Sig,fs,s,stft_type,tf_type,gamma);
STFT_TFD = abs(TFx');
figure; x1=7.5; y1=23; x2=10; y2=33;
set(gcf,'Position',[20 100 640 500]);
imagesc(t,f,STFT_TFD);
axis xy
xlabel('Time/s','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency/Hz','FontSize',20,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'YDir','normal')
%%zoom area
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',1);
h1=axes('position',[0.6 0.2 0.3 0.3]);
axis(h1);
imagesc(t,f,STFT_TFD);
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]); %%不显示刻度
set(gca,'fontsize',12,'linewidth',1,'XColor','r','YColor','r') %%XY坐标轴的颜色
colormap turbo %% 颜色图的色阶

%% TSET1
tf_type ='TSET1';s = 1e-1;
[~,TFx,~,~,Rep,~,q,~,~] = TFM(Sig,fs,s,stft_type,tf_type,gamma);
TSET1_TFD = abs(TFx');
figure
set(gcf,'Position',[20 100 640 500]);
imagesc(t,f,TSET1_TFD);
axis([0 15 0 50]);
xlabel('Time/s','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'YDir','normal')
%%zoom area
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',1);
h1=axes('position',[0.6 0.2 0.3 0.3]);
axis(h1);
imagesc(t,f,TSET1_TFD);
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]); %%不显示刻度
set(gca,'fontsize',12,'linewidth',1,'XColor','r','YColor','r') %%XY坐标轴的颜色
colormap turbo %% 颜色图的色阶

%% SST1
tf_type ='SST1';s = 0.11;
[~,TFx,~,~,Rep,~,q,~,~] = TFM(Sig,fs,s,stft_type,tf_type,gamma);
TSET1_TFD = abs(TFx');
figure
set(gcf,'Position',[20 100 640 500]);
imagesc(t,f,TSET1_TFD);
axis([0 15 0 50]);
xlabel('Time/s','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'YDir','normal')
%%zoom area
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',1);
h1=axes('position',[0.6 0.2 0.3 0.3]);
axis(h1);
imagesc(t,f,TSET1_TFD);
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]); %%不显示刻度
set(gca,'fontsize',12,'linewidth',1,'XColor','r','YColor','r') %%XY坐标轴的颜色
colormap turbo %% 颜色图的色阶

%% TSST1
tf_type ='TSST1';s = 5e-1;
[~,TFx,~,~,Rep,~,q,~,~] = TFM(Sig,fs,s,stft_type,tf_type,gamma);
TSET1_TFD = abs(TFx');
figure
set(gcf,'Position',[20 100 640 500]);
imagesc(t,f,TSET1_TFD);
axis([0 15 0 50]);
xlabel('Time/s','FontSize',24,'FontName','Times New Roman');
ylabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24);
set(gca,'YDir','normal')
%%zoom area
rectangle('Position',[x1 y1 x2-x1 y2-y1],'EdgeColor','r','Linewidth',1);
h1=axes('position',[0.6 0.2 0.3 0.3]);
axis(h1);
imagesc(t,f,TSET1_TFD);
xlim([x1 x2]);ylim([y1 y2]);
set(gca,'xticklabel',[]);set(gca,'yticklabel',[]); %%不显示刻度
set(gca,'fontsize',12,'linewidth',1,'XColor','r','YColor','r') %%XY坐标轴的颜色
colormap turbo %% 颜色图的色阶

% %% Reconstructed amplitudes
% % Amplitude 1
% figure; 
% plot(f,abs(NGDm1),'b','linewidth',1);
% hold on;
% plot(f,abs(NGDMEeD1),'r','linewidth',1); 
% xlabel('Frequency/Hz','FontName','Times New Roman');
% ylabel('Amplitude','FontName','Times New Roman');
% 
% 
% % Amplitude 2
% figure;
% plot(f,abs(NGDm2),'b','linewidth',1);   % time-domain signal modes
% hold on
% plot(f,abs(NGDMEeD2),'r','linewidth',1);
% xlabel('Frequency/Hz','FontSize',24,'FontName','Times New Roman');
% ylabel('Amplitude','FontSize',24,'FontName','Times New Roman');

%% estimated signal
figure
set(gcf,'Position',[20 100 640 500]);	    
set(gcf,'Color','w'); 
plot(t,Sig,'b','linewidth',2);
hold on;
plot(t,(real(iffteASig1)+real(iffteASig2)),'r','linewidth',2);
xlabel('Time (s)','FontSize',24,'FontName','Times New Roman');
ylabel('Amplitude','FontSize',24,'FontName','Times New Roman');
set(gca,'FontSize',24)
set(gca,'linewidth',2);
set(gca,'FontName','Times New Roman')
axis([0 15 -0.2 0.2]);

%% Reconstructed modes
figure
set(gcf,'Position',[20 100 640 500]);	  
subplot(2,1,1)
set(gcf,'Color','w'); 
plot(t,eSig1,'b','linewidth',2)
hold on
plot(t,Sig1-eSig1,'k--','linewidth',2)	% estimated mode
ylabel('C1','FontSize',24,'FontName','Times New Roman');set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);
axis([0 15 -0.2 0.2])
subplot(2,1,2)    
set(gcf,'Color','w'); 
plot(t,eSig2,'b','linewidth',2)
hold on
plot(t,Sig2-eSig2,'k--','linewidth',2)	% estimated mode
ylabel('C2','FontSize',24,'FontName','Times New Roman');
xlabel('Time / Sec','FontSize',24,'FontName','Times New Roman')
set(gca,'YDir','normal')
set(gca,'FontSize',24);
set(gca,'linewidth',2);	
axis([0 15 -0.2 0.2])
%% EQF results
Mode1_EQF_STFT = 20*log10(norm(real(ifftSig1) - real(iffteASig1),2)/norm(real(ifftSig1),2))
Mode2_EQF_STFT = 20*log10(norm(real(ifftSig2) - real(iffteASig2),2)/norm(real(ifftSig2),2))

Mode1_EQF_NGDME = 20*log10(norm(real(ifftSig1) - real(iffteASig1),2)/norm(real(ifftSig1),2))
Mode2_EQF_NGDME = 20*log10(norm(real(ifftSig2) - real(iffteASig2),2)/norm(real(ifftSig2),2))

GD1_EQF_NGDME = 20*log10(norm(eGDest(1,:,end) - gd1,2)/norm(gd1,2))
GD2_EQF_NGDME = 20*log10(norm(eGDest(2,:,end) - gd2,2)/norm(gd2,2))