function [ASpec tbin] = TFspec(GDmulti,IAmulti,band)
%Construct the time-frequency distribution by the estimated group delays
%and the amplitudes
% band = [min max] time band or time range
tnum = 1024; % the number of the time bins
tbin = linspace(band(1),band(2),tnum);
num = size(GDmulti,1); % the number of the signal modes
N = size(GDmulti,2); 
ASpec = zeros(tnum,N);
delta = floor(tnum*0.1e-2);
for kk = 1:num
    temp = zeros(tnum,N);
    for ii = 1:N
        [~,index] = min(abs(tbin - GDmulti(kk,ii)));
        lindex = max(index-delta,1);
        rindex = min(index+delta,tnum);
        temp(lindex:rindex,ii) = IAmulti(kk,ii);
    end
    ASpec = ASpec + temp;
end
ASpec = ASpec';
        
    


