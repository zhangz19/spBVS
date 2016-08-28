function [] = mainProBVS(id, nchain, niter, burn, thin, spatial)
% MCMC setup
verbose = 0;

%============= total iterations, burn-in period and thin (samples every thin iteration) for each MCMC run

nbin = 1e3; % ifverbose, every nbin iterations summarize the acceptance rates for that batch
repl = str2double(num2str(id)); ch0 = repl;
ncase = ceil(repl/nchain);
repl = repl - (ncase-1)*nchain;
fprintf('case = %d, chain = %d:\n', [ncase, repl])  % partition jobs

% data input and initialize parameters
rng('default'); rng(ch0*8);
[x, E] = readDataBVS(ncase, repl, nchain, spatial);

% run MCMC, store the results
xall = cell(1,(niter-burn)/thin); ps = zeros(1,niter);
births = nan(1,niter); deaths = nan(1,niter);
eta_a = 0; batchLen = 50; batchNum = 0; batchTot = niter/batchLen;  eta_rates = zeros(1, batchTot);
sIter = 1;

tic
for iter = 1:niter
    if verbose==1; fprintf('%6d', x.p); if(~mod(iter,20)); fprintf('\n'); end; end
    
    % -------------- update marginal likelihood
    x = getLoglikeBVS(x, E, 0);
    
    % -------------- update random subset inds (model)
    U0 = rand(1); U1 = log(rand(1));
    if (U0 <= 0.5 || x.p == 0) && (x.p < E.pstar) % add another variable
        i0 = randsample(setdiff(1:E.p, x.inds), 1);
        x1 = x; x1.inds = [x.inds, i0]; x1.p = length(x1.inds);
        x1 = getLoglikeBVS(x1, E, 1);
        logratio = x1.logLik - x.logLik - E.alpha;
        if x.p==0; logratio = logratio + log(.5); end %boundary correction
        MHratio = min(logratio, 0); births(iter) = 0;
        if U1 <= MHratio; x = x1; births(iter) = 1; end
    else % delete an existing variable
        i0 = randsample(x.p,1);
        x1 = x; x1.inds(i0) = []; x1.p = length(x1.inds);
        x1 = getLoglikeBVS(x1, E, 1);
        logratio = x1.logLik - x.logLik + E.alpha;
        if x.p == E.pstar; logratio = logratio + log(2); end %boundary correction
        MHratio = min(logratio,0); deaths(iter) = 0;
        if U1 <= MHratio; x = x1; deaths(iter) = 1; end
    end
    
    % -------------- update variation \tau2 and fixed-effects \beta
    Xtild = [ones(E.N,1), E.X(:,x.inds)];
    Sigma = E.M-x.gamma.*E.W; Lo = chol(Sigma, 'lower');
    L = 0.5*sum((Lo'*(x.eta-Xtild*x.beta)).^2) + 0.5*x.invlambda*sum((x.beta).^2);
    x.invtau2 = gamrnd(E.a_tau2+0.5*(E.N+x.p+1), 1/(E.invb_tau2+L));
    
    Sigma = Xtild'*Sigma; Mu = Sigma*x.eta;
    Sigma = Sigma*Xtild + x.invlambda.*eye(x.p+1);
    Lo = chol(Sigma, 'lower');
    x.beta = Lo'\( 1/sqrt(x.invtau2)*randn(size(Mu)) + Lo\Mu );
    mu = Xtild*x.beta;  u = x.eta - mu;
    
    % -------------- update signal-to-noise \Lambda
    x.invlambda = gamrnd(E.a_lambda+0.5*(x.p+1), 1/(E.invb_lambda+0.5*x.invtau2*sum(x.beta.^2)));
    
    if E.spatial == 1
        % -------------- update spatial dependence \gamma
        P = E.loglike0 + E.gammas*(u'*E.W*u)*(0.5*x.invtau2);
        P = exp(P-max(P))/sum(exp(P-max(P)));
        cump = cumsum([0 P(1:(end-1))]);
        U0 = rand(1);
        i0 = sum(U0 > cump);
        x.gamma = E.gammas(1); if(i0>1); x.gamma = E.gammas(i0-1) + E.gap/P(i0)*(U0-cump(i0)); end
    end
    
    % -------------- update latent parameter \eta for Poisson model (E.fixeta==0) using tMALA
    if E.fixeta == 0
        Lo = sqrt(x.invtau2) * chol(E.M - x.gamma*E.W, 'lower');
        v = Lo'*u;
        grad = E.Y - E.offset.*exp(x.eta); if numel(grad>E.H)>0; grad(grad>E.H) = E.H; end
        grad = Lo\grad - v;
        v1 = v + 0.5*E.h*grad + E.h_sqrt*randn([E.N, 1]);
        eta1 = mu + Lo'\v1;
        grad1 = E.Y - E.offset.*exp(eta1); if numel(grad1>E.H)>0; grad1(grad1>E.H) = E.H; end
        grad1 = Lo\grad1 - v1;
        P = sum(E.Y.*eta1 - E.offset.*exp(eta1) - 0.5*v1.^2 + 0.5/E.h*(v1 - v - 0.5*E.h*grad).^2) ...
            - sum(E.Y.*x.eta - E.offset.*exp(x.eta) - 0.5*v.^2 + 0.5/E.h*(v - v1 - 0.5*E.h*grad1).^2);
        U1 = log(rand(1));
        if U1 <= P; x.eta = eta1; eta_a = eta_a + 1; end
    end
    
    if verbose == 1 && iter > nbin && ~mod(iter,nbin)
        rbirths = births((iter-nbin):iter); rdeaths = deaths((iter-nbin):iter);
        fprintf('iter=%d, birth=%.3f, death=%.3f, p=%d, ML = %4.2f, tau2 = %4.2f, lambda = %4.2f\n',...
            [iter, mean(rbirths(~isnan(rbirths))), mean(rdeaths(~isnan(rdeaths))), x.p, x.logLik, 1/x.invtau2, 1/x.invlambda]);
    end
    
    if ~mod(iter, batchLen); batchNum = batchNum+1; eta_a = eta_a/batchLen; eta_rates(batchNum) = eta_a; eta_a = 0; end %disp(eta_rates(batchNum));
    
    % record results
    if(iter > burn) && ~mod(iter, thin);
        xall{sIter} = x; ps(sIter) = x.p;  sIter = sIter+ 1;
    end
end

runtime = toc/60;
fprintf('%d iterations are done with elapsed time %.2f minutes.\n', niter, runtime)
nam = strcat('out',num2str(ncase), '_', num2str(repl),'.mat');
save(nam,'xall','runtime','births','deaths','ps','eta_rates')
end
