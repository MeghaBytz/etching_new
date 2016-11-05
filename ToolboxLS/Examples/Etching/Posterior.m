function [F]  = Posterior(current,index)
global posteriorRecord
[L] = Likelihood(current);
posterior = L + Prior(current,index);
F = posterior;
posteriorRecord = [posteriorRecord posterior];
end
