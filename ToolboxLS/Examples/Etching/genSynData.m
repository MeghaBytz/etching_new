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
%---------------------------------------------------------------------------
% create options structure
[g, options] = createOptionsDefault();

% set case for velocities
options.etchShape = 'varyAlpha';
%---------------------------------------------------------------------------
% Modify the grid.
g.dim = 2;
g.min = [-1; -1];
g.dx = 1 / 100;
g.max = [1;1];
g.bdry = @addGhostExtrapolate;

g = processGrid(g);
%---------------------------------------------------------------------------
% create geometry, layers and mask

% set to zero for vertical stacks
options.horizontal = 0;

% get layer boundaries
% boundaries must be set in getMaterialMap. Haven't figured out an easier
% way yet
[startMatInd, endMatInd, ~, options.map] = getMaterialMap(g, options);
options.colors = getColors(options.map);
options.layer_boundaries = cat(1,startMatInd,endMatInd);
%---------------------------------------------------------------------------
% Plotting parameters
% To plot or not to plot
options.doDisplay = 1;
% Delete previous plot before showing next?
options.deleteLastPlot = 0;
% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
options.useSubplots = 1;

%---------------------------------------------------------------------------
% Integration parameters.
options.tMax = 2;                  % End time.
options.plotSteps = 9;               % How many intermediate plots to produce?
options.t0 = 0;                      % Start time.
options.singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
options.tPlot = (options.tMax - options.t0) / (options.plotSteps - 1);
%---------------------------------------------------------------------------

for i = 1:length(expConditions)
 [data, g, algernonTrenchDim(i,1), algernonTrenchDim(i,2)] = etchingWithOptions(g,options,algernonUnknowns,i);
 kepler(i,:) = algernonTrenchDim(i,:) + abs(normrnd(0,noise,[1 2]));
end

