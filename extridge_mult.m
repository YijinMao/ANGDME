function [tidexmult, tfdv] = extridge_mult(Sig, SampFreq, num, delta, beta,bw,Nfrebin,window)
% Extract ridges for multi-component signals.
% In each iteration,the signal component associated with the extrated ridge is
% reconstructed by a time-frequency (TF) filter and then removed from the original signal so
% that the ridge curves of other signal components with smaller energies
% can be extracted in the subsequent iterations.
%%%%%%%%%%%%%%%%%%%%%%%    input      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sig: FFt of the signal,a row vector
% SampFreq: sampling frequency
% num: the number of the signal components
% delta: maximum allowable time variation between two consecutive points for ridge detection
% beta: controls the smooth degree; the curves will be smoother if beta is smaller
% bw: the time bandwidth of the TF filter (unit£ºsec);
% Nfrebin,window are two parameters for implementing the SFFT


%%%%%%%%%%%%%%%%%%%%%%%%%%%%  output   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fidexmult: obtained ridge index for multi-component signals
% tfdv: the corresponding ridge magnitude 

% if (isreal(Sig))
% Sig = hilbert(Sig);
% end

tidexmult = zeros(num,length(Sig));
tfdv = zeros(num,length(Sig));

for i = 1:num
    [Spec,~,t] = SFFT(Sig(:),SampFreq,Nfrebin,window); % SFFT
    c = findridges(Spec.',delta);
    [sGD,extr_Sig] = Dechirp_filter(Sig,SampFreq,bw,t(c),beta); % the TF filter; we smooth the extracted ridges and the smoothed curve is sGD
    tindex = zeros(1,length(Sig));
    for j = 1:length(Sig)
        [~,tindex(j)] = min(abs(t - sGD(1,j)));
        tfdv(i,j) = abs(Spec(j,tindex(j)));
    end
    tidexmult(i,:) = tindex;
    Sig = Sig - extr_Sig; % remove the extracted signal component so that other ridge curves with smaller energies can be extracted in the subsequent iterations
end

end