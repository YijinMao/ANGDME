function [GDest, Modeest, x] = ANGDME(y, T, eGDset, gamma, tol)

% This code implements the Adptive Nonlinear Group Delay Mode Estimation (ANGDME) algorithm 
%
% This code is based on the GDMD code available from https://www.researchgate.net/publication/345173909_Codes_of_GDMDzip from the following paper
% [1] Chen S, Wang K, Peng Z, et al, Generalized dispersive mode decomposition: Algorithm and applications, Journal of Sound and Vibration, 2020.
% 
% Please check the accompanying license and the license of [1] before using. 
%
%
% Inputs:
%    y:  unilateral Fourier spectrum of the measured signal, a row vector 
%    T:  time duration 
%    eGDset:  initial group delay for the signal modes, a matrix (each row corresponds to the initial group delay of a certain signal mode)
%    lambda:  parameter controlling the sparsity --> adptively
%    gamma:   parameter controling the smooth degree of group delay increment
%    tol:     tolerance of convergence criterion
% Outputs:
%    GDest:  the estimated group delays at each iteration 
%    sest:   the estimated signal modes in the frequency domain at each iteration 
%
%
% ANGDME Authors: Yijin Mao, Hao Liang (haoliang@stu.xmu.edu.cn) and Xiaotong Tu (xttu@xmu.edu.cn)
% Last modified by: 22/08/20

% parameter setting
[K, N] = size(eGDset); % N is the number of the frequency samples, K is the number of the signal components; 
f = (0:N-1)/T; 
e = ones(N,1); e2 = -2*e;
H = spdiags([e e2 e], 0:2, N-2, N); % the second-order difference matrix
HtH = H'*H; 
tempm = repmat('H,', 1, K);
D = eval(sprintf('blkdiag(%s)',tempm(1:end-1)));
A = spdiags((zeros(N,1)), 0, N, N*K);

iternum = 500;  % the maximum allowable iterations
GDiter = zeros(K, N, iternum); 
siter = zeros(K, N, iternum); 

% start iteration
iter = 1; sDif = tol + 1;
while ( sDif > tol && iter <= iternum ) 
    for kk = 1:K 
        expm = exp(-1i*(2*pi*(cumtrapz(f,eGDset(kk,:)))));
        Aerm = spdiags(expm(:), 0, N, N);
        A(:,((kk-1)*N+1):N*kk) = Aerm;  % update A
    end
% update demodulated signals 
     x = Estimate_NCS(A,D,y(:),K,N); % update X

    for kk = 1:K
        xm = x(((kk-1)*N+1):N*kk);
        deltaphase = unwrap(angle(xm));  % get the phase
        deltaGD = Differ(deltaphase,1/T)/2/pi;
        deltaGD = (gamma*HtH + speye(N))\deltaGD(:);
        eGDset(kk,:) = abs(eGDset(kk,:) - deltaGD.'); % update the group delay \tao k
        siter(kk,:,iter) =  A(:,((kk-1)*N+1):N*kk)*xm; % update Yk
    end
    GDiter(:,:,iter) = eGDset;
    
    if iter>1
        sDif = 0;
        for kk = 1:K
            sDif = sDif + (norm(siter(kk,:,iter) - siter(kk,:,iter-1))/norm(siter(kk,:,iter-1))).^2;
        end
    end
    iter
    sDif
    
    if iter<20
        sDif = 1;
    end
    iter = iter + 1;
end
    GDest = GDiter(:,:,1:iter-1);     % group delay {\tao k}
    Modeest = siter(:,:,1:iter-1);    % signal components {Yk}
end




    
    
    
    