clear all
global expConditions
%global data
global noUnknowns
global numberOfIons
global numberOfRadicals
global noExperimentalMetrics
global levelSetUnknowns
global numberOfMaterials

noExperimentalMetrics = 2;
numberOfIons = 2;
numberOfRadicals = 2;
levelSetUnknowns = 2;
numberOfMaterials = 1;
noUnknowns = numberOfIons*2+numberOfRadicals*2+levelSetUnknowns+1;

load allSynExpConditions

expConditions = allSynExpConditions;

%make fake real parameters for novel resist material "Algernon"
A1 = [.9 2];
B1 = [3.1 5.7];
rxnProb1 = [0.4 .7];
SCl1 = [.4 .7];
lambda1 = 1;
lambda2 = .6;
epsS1 = 50;
epsD1 = 8;
noise = .02;

A2 = [.2 .1];
B2 = [1 1];
rxnProb2 = [0.4 .8];
SCl2 = [.9 .9];
epsS2 = 11;
epsD2 = 2;

material1 = [A1 B1 rxnProb1 SCl1 lambda1 lambda2 epsS1 epsD1 noise];
material2 = [A2 B2 rxnProb2 SCl2 lambda1 lambda2 epsS2 epsD2 noise];

algernonUnknowns = [material1; material2]

%data = allSynData(y(1:noTrainingData));
for i = 6%length(expConditions)
 [data, g, algernonTrenchDim(i,1), algernonTrenchDim(i,2)] = etchingWithMask('low', 'contour', 'varyAlpha',algernonUnknowns,i);
 kepler(i,:) = algernonTrenchDim(i,:) + abs(normrnd(0,noise,[1 2]));
end

