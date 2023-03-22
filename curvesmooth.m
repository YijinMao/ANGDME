function outf = curvesmooth(f,beta)
% f: input instantaneous frequencies (IFs) or group delays (GDs)
% outf: output smoothed IF or GD curves
% beta: controls the smooth degree; the curves will be smoother if beta is smaller

[K,N] = size(f); 
e = ones(N,1); e2 = -2*e; % e2(1) = -1;e2(end) = -1;
oper = spdiags([e e2 e], 0:2, N-2, N);
opedoub = oper'*oper;
outf = zeros (K,N);
for i = 1:K
    outf(i,:) = (2/beta*opedoub + speye(N))\f(i,:).'; % smooth the instantaneous group delay curves by low pass filtering
end

end