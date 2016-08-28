function [x, E] = readDataBVS(ncase, repl, nchain, spatial)
% import data, initialize parameter, set environmental variables
load('yourdata.mat') %X full design matrix, Y response, W spatial adjacency
E.N = size(X, 1); %size(W,1);
E.X = zscore(X);  %standardized predictors without intercept
E.Y = Y;
E.p = size(E.X, 2);
E.pstar = min(E.p, E.N-1); %upper bound of p
E.alpha = 0; % Penalty term: 1: AIC, log(n)/2: BIC
E.spatial = spatial;  % spatial or nonspatial?
if E.spatial == 1
    E.W = W;   E.M = diag(sum(W,1)); invM = inv(E.M);
    eigs = eig(sqrt(invM)*W*sqrt(invM));
    lgamma = max(1/min(eigs),-1); ugamma = 1/max(eigs);
    gap = 1e-3; gammas = (lgamma+gap):gap:(ugamma-gap); len = length(gammas);
    if(~exist('loglike0.mat','file')) % precalculaiton: make CAR model sampling much faster!
        loglike0 = zeros(1,len); for i = 1:len; loglike0(i) = 0.5*sum(log(eig(E.M-gammas(i)*E.W))); end; save('loglike0.mat','loglike0');
    end
    load('loglike0.mat')
    E.gap = gap; E.gammas = gammas; E.loglike0 = loglike0;
else
    E.W = 0; E.M = eye(E.N);
end

% Poisson part
E.fixeta = 1; % 1 = Gaussian continuous response, 0 = Poisson count response
if ncase > 2; E.fixeta = 0; end

if E.fixeta == 0
    E.offset =  ones(E.N,1);     %no offset
    % % for truncated MALA
    E.H = 50; E.h = 0.0016; E.h_sqrt = sqrt(E.h);
end

% set prior
meantau2 = 0.01; vartau2 = 10^2;
E.a_tau2 = 2+meantau2^2/vartau2;
E.invb_tau2 = meantau2*(E.a_tau2-1);

meanlambda = 1; varlambda = 1e2; %100
E.a_lambda = 2+meanlambda^2/varlambda;
E.invb_lambda = meanlambda*(E.a_lambda-1);

% may consider certain transformation here
logTransformY = 0; if logTransformY == 1; x.eta = log(E.Y); else x.eta=E.Y; end
if E.fixeta == 0;  x.eta = log((E.Y+ 0.5*(E.Y==0))./E.offset); end % for Poisson

% set initial value
p0s = floor(linspace(1,E.p,nchain)); % initial number of variable spreads out its support
p0 = p0s(repl);
x.inds = randsample(1:E.p, p0); x.p = length(x.inds);
Xtild = [ones(E.N,1), E.X(:,x.inds)];
x.beta = (Xtild'*Xtild)\(Xtild'*x.eta) .* (1 + 1/50*randn([size(Xtild,2),1]));
x.invtau2 = 1/( sum((sqrt(diag(E.M)).*(x.eta-Xtild*x.beta)).^2)/(E.N-E.p) );

x.gamma = 0;  x.invlambda = 1/meanlambda;
x = getLoglikeBVS(x, E, 1);
end
