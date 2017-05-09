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
numberOfExperiments = length(expConditions(:,1));
%make fake real parameters for novel resist material "Algernon"
A1 = [.1 .2];
B1 = [.13 .05];
rxnProb1 = [0.82 .97];
SCl1 = [.92 .97]; %these were set
epsS1 = 11;
epsD1 = 4;%these were set
lambda1 = 2;
lambda2 = .7;
noise = .02;

material1 = [A1 B1 rxnProb1 SCl1 epsS1 epsD1 lambda1 lambda2 noise];

algernonUnknowns = [material1]
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
options.horizontal = 1;

% get layer boundaries
% boundaries must be set in getMaterialMap. Haven't figured out an easier
% way yet
[startMatInd, endMatInd, ~, options.map] = getMaterialMap(g, options);
options.colors = getColors(options.map);
options.layer_boundaries = cat(1,startMatInd,endMatInd);
%---------------------------------------------------------------------------
% Plotting parameters
% To plot or not to plot
options.doDisplay = 0;
% Delete previous plot before showing next?
options.deleteLastPlot = 0;
% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
options.useSubplots = 1;

%---------------------------------------------------------------------------
% Integration parameters.
options.tMax = 1; % End time.
options.plotSteps = 4;  % How many intermediate plots to produce?
options.t0 = 0; % Start time.
options.singleStep = 0; % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
options.tPlot = (options.tMax - options.t0) / (options.plotSteps - 1);
%---------------------------------------------------------------------------
%etchingWithOptions(g,options,algernonUnknowns,1)
%load thetaSequenceRODEoFinal
%theta = thetaSequenceRODEoFinal(:,:,6);
pd = 10;


for qq=1:numberOfExperiments
    [VxAlgernon, VyAlgernon] = plasma(algernonUnknowns,qq);
end
elapsedtime = toc;


