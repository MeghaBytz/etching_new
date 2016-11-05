function [alpha, t, a, prob,PosteriorCatch] = MetropolisHastings(theta,current,PosteriorCurrent,index)
%independent sampling from proposal function (current is put here to hold
%unused blocks of parameters)
global alphaRecord
[new,chainSD] = ProposalFunction(theta,current,index);
[PosteriorNew] = Posterior(new,index);
[alpha] = exp(PosteriorNew + ProposalPdf(current,new,index,chainSD)-(PosteriorCurrent+ProposalPdf(new,current,index,chainSD))); % Ratio of the density at the candidate (theta_ast) and current (current) points
alphaRecord = [alphaRecord alpha];
if rand <= min(alpha,1)
   t    = new;        % Accept the candidate
   prob = min(alpha,1);     % Accept with probability min(alpha,1)
   a    = 1;                % Note the acceptance
   PosteriorCatch=PosteriorNew;
else
   a = 0;
   t    = current;            % Reject the candidate and use the same state
   prob = 1-min(alpha,1);   % The same state with probability 1-min(alpha,1)
   PosteriorCatch = PosteriorCurrent;
end
end
