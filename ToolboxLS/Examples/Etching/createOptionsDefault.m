function [g, options] = createOptionsDefault()

options.doMask = 0;
options.accuracy = 'medium';
options.etchShape = 'rectangle';
%---------------------------------------------------------------------------
% Create the grid.
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
[startMatInd, endMatInd, options.map] = getMaterialMap(g,options);
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
% Pause after plot?
options.pauseAfterPlot = 0;
%---------------------------------------------------------------------------
% Integration parameters.
options.tMax = 2;                  % End time.
options.plotSteps = 9;               % How many intermediate plots to produce?
options.t0 = 0;                      % Start time.
options.singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
options.tPlot = (options.tMax - options.t0) / (options.plotSteps - 1);
%---------------------------------------------------------------------------