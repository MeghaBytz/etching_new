function [F, chainSD] = ProposalFunction(theta,current,index)
global proposedParameterRecord
global proposalLB
global proposalUB
global proposalMean
global proposalSD
parameter = current;
parameter(index) = (proposalUB(index)-proposalLB(index)).*rand(1,1) + proposalLB(index);
chainSD = 1; %placeholder
F = parameter;
end
