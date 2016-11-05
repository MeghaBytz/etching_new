% syms u v
% 
% S = solve([k1*A*M-k2*A1*M-k3*A1 == 0, u - v == 1], [u, v])

 %Make reaction matrix
 numberOfReactions = 2;
 numberOfSpecies = 5;
 M = zeros(numberOfReactions, numberOfSpecies);
 kRates = [1 1];
 speciesConcentrations = ones(1,numberOfSpecies);
 
 for i = 1:numberOfSpecies
     for j = 1:numberOfReactions
        M(:,i)*kRates
     end
 end
 
 