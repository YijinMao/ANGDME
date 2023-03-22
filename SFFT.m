function [Spec,f,t] = SFFT(sigfft,SampFreq,N,WinLen)

% short-frequency Fourier transform

if(nargin < 3),
    error('At least 3 inputs are required!');
end

% SigLen =N;
% SigLen = 2*length(sigfft)+1;
SigLen = 2*(length(sigfft)-1)+1;
% N = length(sigfft);
if (nargin < 3)
    WinLen = N / 4;
end


%RatioNum = length(Ratio);

%% compute frequency spectrum
%Spec = fft(Sig)./length(Sig); %take the fft of our sin wave, y(t)
%Freq=( -ceil((N-1)/2):N-1-ceil((N-1)/2) )/N*SampFreq;
%df = Freq(ceil(end/2)+1:ceil(end));
N = length(sigfft);
df = [0:(N-1)]/N*SampFreq/2;
%df=[0:(N/2)]*SampFreq/N;
%sigfft = Spec(1:floor(end/2));
%lllen=length(sigfft)

%%
WinLen = ceil(WinLen / 2) * 2;
f = linspace(-1,1,WinLen)';
WinFun = exp(log(0.005) * f.^2 );
WinFun = WinFun / norm(WinFun);
Lh = (WinLen - 1)/2; 

fft_len = length(df);
Spec = zeros(SigLen,fft_len) ;   % matrix
conjSpec = Spec;

%%

%wait = waitbar(0,'Please wait...');
for iLoop = 1:fft_len
    
    %waitbar(iLoop/fft_len,wait);
   
    tau = -min([round(fft_len/2)-1,Lh,iLoop-1]):min([round(fft_len/2)-1,Lh,fft_len-iLoop]);  % signal span
    temp = floor(iLoop + tau);
    sSig = sigfft(temp);

    temp1 = floor(Lh+1+tau);    % window span
    sSig = sSig .* conj(WinFun(temp1)); % Z(t)* complex conjugate of window?
    Spec(1:length(sSig),iLoop) = sSig;  % windowed analytic signal
    conjSpec(:,iLoop)  =  fliplr(Spec(:,iLoop));
end
%%
iSpec  = [conjSpec(1:end-1,:);Spec];
iLen = length(iSpec);
Spec = iSpec(round(iLen/2):end,:);

Spec = ifft(Spec);
Spec = abs(Spec)/2/pi;

%close(wait);

[SigLen,nLevel] = size(Spec);

f = [0:nLevel-1]/nLevel * SampFreq/2;  % frequency  in TF plane?
t = (0: SigLen-1)/SampFreq;      % time  in TF plane
Spec =Spec';
%[fmax fmin] = FreqRange(Sig);
%fmax = fmax * SampFreq;
%fmin = fmin * SampFreq;
% mesh(t,f,Spec);  
% % axis([min(t) max(t) fmin fmax]);
% xlabel('Time / us','FontName','Arial','FontSize',26)  
% ylabel('Freq / MHz ','FontName','Arial','FontSize',26)
% set(gca,'FontName','Arial','FontSize',26)
% set(gcf,'position',[30 30 640 500])
% colormap('jet')
%=====================================================%
% Plot the result                                     %
%=====================================================%

end

