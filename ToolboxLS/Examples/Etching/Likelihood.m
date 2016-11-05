function [F] = Likelihood(current)
global expData;
global likelihoodRecord
global etchRecord
global noExperimentalMetrics
likeData = 0;
noise = .02;%current(end); %change this back to unknown noise parameter eventually
trenchDim = zeros(noExperimentalMetrics,1);
for i=1:length(expData)
  [a, g, trenchDim(1), trenchDim(2)] = etchingWithMask('low', 'contour', 'varyAlpha',current,i)
  trenchDim = abs(normrnd(0,current(end),[2 1]))+trenchDim;
   for j = 1:length(noExperimentalMetrics)
            Like = normpdf(expData(i,j),trenchDim(j),noise) + 10e-20; %change back to look at both horizontal and vertical etch rates according to data function
            likeData = log(Like) + likeData;
    end
end
likelihoodRecord = [likelihoodRecord likeData];
F = likeData;
