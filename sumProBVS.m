%% Summarize Bayesian variable selection output
function [] = sumProBVSs(spatial)
% summary function for BVS results
collectit = 1;
for ncase = 1:1
    % fprintf('case = %d\n', ncase)
    if collectit == 1
        load('yourdata.mat')
        E.X = zscore(X);  %standardized predictors without intercept
        [~, p] = size(E.X);
        E.spatial = 0;
        
        dirs = dir(strcat('out',num2str(ncase),'_*'));
        chs = 1:length(dirs); nch = length(chs);
        for ch = 1:nch
            load(dirs(chs(ch)).name)
            if ch == 1 % initialization
                niter = length(xall); burn = 0; thin = 1; realburn = 0; %#ok<*USENS>
                nsample = (niter-burn)/thin;    tot = nch*nsample;
                allps = zeros(nsample,nch); allgammas = zeros(nsample,nch); alltau2s = zeros(nsample,nch); alllambdas = zeros(nsample,nch);
                allbirths = zeros(nsample,nch); alldeaths = zeros(nsample,nch);
                indmat = zeros(tot,p); betamat = zeros(tot,p); beta0mat = zeros(tot,1);
            end
            allps(:,ch) = ps(realburn+burn+(1:nsample).*thin);
            allbirths(:,ch) = births(realburn+burn+(1:nsample).*thin);
            alldeaths(:,ch) = deaths(realburn+burn+(1:nsample).*thin);
            for i = 1:nsample
                indmat((ch-1)*nsample+i,xall{burn+i*thin}.inds) = 1;
                beta0mat((ch-1)*nsample+i) = xall{burn+i*thin}.beta(1);
                betamat((ch-1)*nsample+i,xall{burn+i*thin}.inds) = xall{burn+i*thin}.beta(2:end)';
                allgammas(i,ch) = xall{burn+i*thin}.gamma;
                alltau2s(i,ch) = 1/xall{burn+i*thin}.invtau2;
                alllambdas(i,ch) = 1/xall{burn+i*thin}.invlambda;
            end
        end
        save(strcat('sumout',num2str(ncase),'.mat'),'allps','tot','allbirths','alldeaths','indmat','betamat','allgammas','alltau2s')
    end
    load(strcat('sumout',num2str(ncase),'.mat'))
    kplot = 1;
    subplot(2,2,kplot); plot(allps); kplot = kplot+1;
    title('Traceplot of the number of selected variables','FontSize',8)
    ps = reshape(allps, [1,tot]);
    fprintf('\nPosterior probability mass of the number of selected variables:\n')
    tab = tabulate(ps); 
    fprintf('Value    Count    Percent\n')
    disp(num2str(tab(tab(:,2)>0,:)))
    mbirth = nanmean(reshape(allbirths, [1,tot]));
    mdeath = nanmean(reshape(alldeaths, [1,tot]));
    fprintf('\nmean birth rate = %5.4f, death rate = %5.4f.\n\n', [mbirth, mdeath])
    
    tau2s = reshape(alltau2s, [1,tot]);
    subplot(2,2,kplot); plot(alltau2s); kplot = kplot+1;
    title('Traceplot of the noise level','FontSize',8)
    disp('Posterior mean [lower, upper] 95% credible interval for noise level:')
    %[l,u] = FindHPDset(tau2s', 0.95, []);
    l = quantile(tau2s',0.025);
    u = quantile(tau2s',0.975);
    fprintf('  %5.4f  [%5.4f, %5.4f]\n\n', [mean(tau2s),l,u]); 
    
    lambdas = reshape(alllambdas, [1,tot]);
    subplot(2,2,kplot); plot(alllambdas); kplot = kplot+1;
    title('Traceplot of the overall sigma-to-noise ratio','FontSize',8)
    disp('Posterior mean [lower, upper] 95% credible interval for overall sigma-to-noise ratio:')
    %[l,u] = FindHPDset(tau2s', 0.95, []);
    l = quantile(lambdas',0.025);
    u = quantile(lambdas',0.975);
    fprintf('  %5.4f  [%5.4f, %5.4f]\n\n', [mean(lambdas),l,u]); 
    
    if spatial == 1
        gammas = reshape(allgammas, [1,tot]);
        subplot(2,2,kplot); plot(allgammas);
        title('traceplot of the spatial dependence','FontSize',8)
        disp('Posterior mean [lower, upper] 95% credible interval for spatial dependence:')
        %[l,u] = FindHPDset(gammas', 0.95, []);
        l = quantile(gammas',0.025);
        u = quantile(gammas',0.975);
        %disp([mean(gammas),l,u]);
        fprintf('  %5.4f  [%5.4f, %5.4f]\n\n', [mean(gammas),l,u]); 
    end
    
    % -------------- variable wise summary
    P = size(betamat, 2);
    mat = zeros(P, 3);
    for i = 1:P
        betavec = betamat(:,i); betavec = betavec(betavec~=0);
        mat(i,:) = [mean(betavec), quantile(betavec,0.025), quantile(betavec, 0.975)];
    end
    mat_var = mean(indmat); np = min(15, length(mat_var));
    disp('Top variable index (marginal inclusion probability), posterior mean [lower, upper] 95% CI of effect')
    [mat_var,I] = sort(mat_var,'descend'); mat_var = [I',mat(I,:), mat_var']; %disp(num2str(mat_var(1:np,:))); fprintf('\n')
    for i = 1:np; fprintf('%2d (%5.4f):   %5.4f  [%5.4f, %5.4f]\n', mat_var(i,[1,5,2:4]));  end; fprintf('\n')
    
    % -------------- model wise summary
    [xu, ~, k] = unique(indmat, 'rows');
    count = histc(k, 1:size(xu, 1));
    [~,I] = sort(count,'descend');
    tab = [xu, count/sum(count)];
    tab = tab(I,:);
    bmat = nan(size(tab,1),size(tab,2)-1);
    lb = bmat; ub = bmat;
    bmat0 = nan(size(tab,1),1); lb0 = bmat0; ub0 = bmat0;
    for i = 1:length(I)
        bmat(i,:) = mean(betamat(k==I(i), :),1);
        lb(i,:) = quantile(betamat(k==I(i), :),0.025,1);
        ub(i,:) = quantile(betamat(k==I(i), :),0.975,1);
        bmat0(i,:) = mean(beta0mat(k==I(i)));
        lb0(i,:) = quantile(beta0mat(k==I(i)),0.025);
        ub0(i,:) = quantile(beta0mat(k==I(i)),0.975);
    end
    save(strcat('paras',num2str(ncase),'.mat'),'allps', 'ps','tab','bmat','lb','ub','mat_var','bmat0','lb0','ub0')
end
end
