function [sGD, out_Sig] = Dechirp_filter(Sig, SampFreq, bw, GD, beta)
% Extract a signal component by using de-chirp and low-pass filtering techniques
% This technique can be regarded as a time-frequency (TF) filtering technique
% Sig: fft of the signal, a row vector
% SampFreq: sampling frequency
% bw: bandwidth of the TF filter in time domain (unit sec)
% GD: group delay series of the target component£»it can be extracted by detecting TF ridges
% beta: controls the smooth degree the GD; the curves will be smoother if beta is smaller
% out_Sig: extracted signal in frequency domain
% sGD: smoothed GD

% if (isreal(Sig))
% Sig = hilbert(Sig);
% end

Nf = length(Sig);
f = [0:Nf-1]/Nf * SampFreq/2;
sGD = curvesmooth(GD,beta); % smooth the GD
phase = cumtrapz(f,sGD); % integral
Sig1r = Sig .* exp( 1j * 2 * pi * phase);   % de-chirp

T = 2*Nf/SampFreq; % time duration
Sig_comp = low_filter(Sig1r,T,bw/2); % low pass filtering
out_Sig = Sig_comp.* exp(- 1j * 2 * pi * phase); % recover the filtered signal

end
