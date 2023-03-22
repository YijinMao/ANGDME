function [weights,lambda_est,ML] = MCMCS_LSM(PHI,t,a,b,lambda_init,eta,cd)
%---------------------------------------------------------------------------------------
% This Algorithm is for the following manuscript
% Complex Multitask Bayesian Compressive Sensing Using Laplacian Scale Mixture Prior
% authored by Qilei Zhang & Yu Lei Email:zhangqilei@nudt.edu.cn
% The Code is based on the MT-CS code and FastLaplace code available on the
% corresponding sites
% Modified by Yu Lei in 2021.10.25
%---------------------------------------------------------------------------------------
% The MT-CS algorithm for the following paper:
% "Multi-Task Compressive Sesning" (Preprint, 2007). The algorithm 
% is an extension of the fast RVM algorithm [Tipping & Faul, 2003]
% in two-fold: (i) the noise variance is marginalized and (ii) it is for
% multi-task CS, including single-task CS as a special case
% Coded by: Shihao Ji, ECE, Duke University
% last change: May. 15, 2007
% A bug was fixed on Aug.03, 2008 for the cases where signals are dramatic undersampled
%---------------------------------------------------------------------------------------
% The implements the fast Laplace algorithm from the following paper:
% [1] S. D. Babacan, R. Molina, A. K. Katsaggelos. Bayesian Compressive Sensing using Laplace Priors,
% submitted for publication, IEEE Transactions on Image Processing, September 2008.
%% ---------------------------------------------------------------------------------------
% Input:
%   PHI: projection matrix. Cell structure, One cell for one task.
%   t:   CS measurements. Cell structure, One cell for one task.
%   a,b: parameters of Gamma prior on noise variance
%   lambda_init: initial value of lambda. [] is recommended and means that the value of lambda is autoestimated.
%   eta: threshold for stopping the algorithm (suggested value: 1e-8)
% Output:
%   weights: sparse weights for all the tasks. One column for one task
%   lambda_est: estimated value of lambda
%   ML:      the increase of the joint mariginal likelihood for each iteration
%%
%----------------------------------------------
% Initialized value of c & d
if nargin < 7    
    cd=[1 0.5];
end
c=cd(1);d=cd(2);
if nargin < 5,
    eta = 1e-8;
end
if nargin < 4,
    lambda_init = [];
end
%
% Obtain the task number 
if iscell(t)
    NT = length(t);
else
    NT = 1;
    PHI = {PHI};
    t = {t};
end
%fprintf(1,'This is a %d-task learning!\n',NT);
    
% find initial alpha
for k = 1:NT
    [N(k),M(k)] = size(PHI{k});
end
if sum(abs(M-M(1))) ~= 0
    error('Sorry! The sizes of the underlying signals should be the same!\n');
else
    M = M(1);
end
% find initial alpha
K = repmat(N+a,[M,1]);
%
for k = 1:NT
    PHIt(:,k) = PHI{k}'*t{k};
    PHI2(:,k) = sum(PHI{k}.*conj(PHI{k}))';
    G2(k) = t{k}'*t{k}+b;
end
G2 = repmat(G2,[M,1]);
X = G2.*PHI2./(PHIt.*conj(PHIt));
ml = K.*log(X./K)-(K-1).*log((X-1)./(K-1));

ml_sum = sum(ml,2);
while 1
    [ML,index] = max(ml_sum);
    alpha = NT./real(sum((K(index,:).*(PHIt(index,:).*conj(PHIt(index,:)))./G2(index,:)-PHI2(index,:))./...
        (PHI2(index,:).*(PHI2(index,:)-(PHIt(index,:).*conj(PHIt(index,:)))./G2(index,:))),2));
    if alpha > 0  % alpha should be greater than 0
        break;
    else
        ml_sum(index) = 0;
    end
end

for k = 1:NT
    % compute initial mu, Sig, S, Q, G
    phi{k} = PHI{k}(:,index);
    Hessian = real(alpha+phi{k}'*phi{k});
    Sig{k} = 1/Hessian;
    mu{k} = Sig{k}*PHIt(index,k);
    left = PHI{k}'*phi{k};
    S(:,k) = real(PHI2(:,k)-Sig{k}*left.*conj(left));
    Q(:,k) = PHIt(:,k)-Sig{k}*PHIt(index,k)*left;
    G(:,k) = real(G2(:,k)-Sig{k}*(PHIt(index,k).*conj(PHIt(index,k))));
end
clear PHI2 left;
%
selected = index;
deleted = [];

for count = 2:10000
    
    %index
    
    s = S; q = Q; g = G;
    Alpha = repmat(alpha,[1,NT]);
    s(index,:) = Alpha.*S(index,:)./(Alpha-S(index,:));
    q(index,:) = Alpha.*Q(index,:)./(Alpha-S(index,:));
    g(index,:) = g(index,:)+Q(index,:).*conj(Q(index,:))./(Alpha-S(index,:));
    
    Delta_M=sum((K.*(q.*conj(q))./g-s)./(s.*(s-(q.*conj(q))./g)),2);
    %theta = NT./sum((K.*(q.*conj(q))./g-s)./(s.*(s-(q.*conj(q))./g)),2);
    
    if isempty(lambda_init),
        lambda = ( length(index) - 1 + c ) / (sum(1./alpha) + d); % lambda in (24)
    else
        lambda = lambda_init;
    end
    
    NewAlphas=(NT+sqrt(NT^2+4*lambda.*Delta_M))./Delta_M/2;
    
    % choice the next alpha that maximizes marginal likelihood
    ml  = repmat(-inf,[M,1]);
    ig0 = find(Delta_M>0);
    % index for re-estimate
    [ire,~,which] = intersect(ig0,index);
    if ~isempty(ire)
        Alpha1 = repmat(NewAlphas(ire),[1,NT]);
        Alpha0 = repmat(alpha(which),[1,NT]);
%         delta = 1./Alpha1-1./Alpha0;
        
        ml_New=-lambda./NewAlphas(ire)+sum(log(Alpha1./(Alpha1+s(ire,:))),2)...
            -sum(K(ire,:).*log(1-(q(ire,:).*conj(q(ire,:)))./(g(ire,:).*(Alpha1+s(ire,:)))),2);
        ml_Old=-lambda./alpha(which)+sum(log(Alpha0./(Alpha0+s(ire,:))),2)...
            -sum(K(ire,:).*log(1-(q(ire,:).*conj(q(ire,:)))./(g(ire,:).*(Alpha0+s(ire,:)))),2);
        ml(ire)=ml_New-ml_Old;
        
    end
    % index for adding
    iad = setdiff(ig0,ire);
    if ~isempty(iad)
        Alpha2 = repmat(NewAlphas(iad),[1,NT]);
        ml(iad)=-lambda./NewAlphas(iad)+sum(log(Alpha2./(Alpha2+s(iad,:))),2)...
            -sum(K(iad,:).*log(1-(q(iad,:).*conj(q(iad,:)))./(g(iad,:).*(Alpha2+s(iad,:)))),2);
        which = intersect(deleted,iad);
        ml(which) = -inf;
    end
    is0 = setdiff([1:M],ig0);
    % index for deleting
    [ide,~,which] = intersect(is0,index);
    if ~isempty(ide)
         Alpha3 = repmat(alpha(which),[1,NT]);
         if length(index) == 1,
             ml(ide) = -inf;
         else
             ml(ide) =lambda./alpha(which)-sum(log(Alpha3./(Alpha3+s(ide,:))),2)...
            +sum(K(ide,:).*log(1-(q(ide,:).*conj(q(ide,:)))./(g(ide,:).*(Alpha3+s(ide,:)))),2);
         end
               
    end
    [ML(count),idx] = max(real(ml));

    % check if terminates?
    if count > 2 & abs(ML(count)-ML(count-1)) < (max(ML)-ML(count))*eta
        break;
    end
    % update alphas
    which = find(index==idx);
    if Delta_M(idx) > 0
        if ~isempty(which) % re-estimate
            Alpha = NewAlphas(idx);
            delta = Alpha-alpha(which);
            for k = 1:NT
                Sigii = Sig{k}(which,which); mui = mu{k}(which); Sigi = Sig{k}(:,which);
                ki = delta/(1+Sigii*delta);
                mu{k} = mu{k}-ki*mui*Sigi;
                Sig{k} = Sig{k}-ki*Sigi*Sigi';
                comm = PHI{k}'*(phi{k}*Sigi);
                S(:,k) = real(S(:,k) + ki*(comm.*conj(comm)));
                Q(:,k) = Q(:,k) + ki*mui*comm;
                G(:,k) = real(G(:,k) + ki*(Sigi'*PHIt(index,k))*conj(Sigi'*PHIt(index,k)));
            end
            %
            alpha(which) = Alpha;
        else % adding
            Alpha = NewAlphas(idx);
            for k = 1:NT
                phii = PHI{k}(:,idx); Sigii = 1/(Alpha+S(idx,k)); mui = Sigii*Q(idx,k);
                comm1 = Sig{k}*(phi{k}'*phii);
                ei = phii-phi{k}*comm1;
                off = -Sigii*comm1;
                Sig{k} = [Sig{k}+Sigii*comm1*comm1', off; off', Sigii];
                mu{k} = [mu{k}-mui*comm1; mui];
                comm2 = PHI{k}'*ei;
                S(:,k) = real(S(:,k) - Sigii*(comm2.*conj(comm2)));
                Q(:,k) = Q(:,k) - mui*comm2;
                G(:,k) = real(G(:,k) - Sigii*(t{k}'*ei).*conj(t{k}'*ei));
                phi{k} = [phi{k},phii];
            end
            %
            index = [index;idx];
            alpha = [alpha;Alpha];
        end
    else
        if ~isempty(which) % deleting
            for k = 1:NT
                Sigii = Sig{k}(which,which); mui = mu{k}(which); Sigi = Sig{k}(:,which);
                Sig{k} = Sig{k}-Sigi*Sigi'/Sigii; Sig{k}(:,which) = []; Sig{k}(which,:) = [];
                mu{k}  = mu{k}-mui/Sigii*Sigi; mu{k}(which) = [];
                comm = PHI{k}'*(phi{k}*Sigi);
                S(:,k) = real(S(:,k) + (comm.*conj(comm))/Sigii);
                Q(:,k) = Q(:,k) + mui/Sigii*comm;
                G(:,k) = real(G(:,k) + (Sigi'*PHIt(index,k)).*conj(Sigi'*PHIt(index,k))/Sigii);
                phi{k}(:,which) = [];
            end
            %
            index(which) = [];
            alpha(which) = [];
        end
    end

end
% output
weights	= zeros(M,NT);
for k = 1:NT
    weights(index,k) = mu{k};
end
lambda_est=lambda;



