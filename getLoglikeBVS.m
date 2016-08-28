function [x] = getLoglikeBVS(x, E, proposeParameter)
Xtild = [ones(E.N,1), E.X(:, x.inds)]; % include intercept
Lo = E.M-x.gamma.*E.W; L = x.eta'*Lo*x.eta;
Lo = Xtild'*Lo; Mu = Lo*x.eta;
Lo = Lo*Xtild + x.invlambda.*eye(x.p+1);
Lo = chol(Lo, 'lower');
Mu = Lo\Mu; L = (L - sum(Mu.^2))/2;
n1 = 0.5*E.N;
if proposeParameter == 1
    x.invtau2 = gamrnd(E.a_tau2+n1, 1/(E.invb_tau2+L));
    Mu = Mu + sqrt(x.invtau2)*randn(size(Mu));
    x.beta = Lo'\Mu;
end
x.logLik = 0.5*(x.p+1)*log(x.invlambda) - sum(log(diag(Lo))) - (E.a_tau2+n1)*log(1+L/E.invb_tau2);
end
