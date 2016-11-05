function [ data, g, data0 ] = etchingVelocityFunc(accuracy, displayType, etchShape)
% etchingTest: demonstrate simple etching, using a constant etching rate.
%
%   [ data, g, data0 ] = etchingTest(accuracy, displayType)
%  
%
% Parameters:
%
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst (default).
%                  'medium'      Use odeCFL2 and upwindFirstENO2.
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   displayType  String to specify how to display results.
%                  The specific string depends on the grid dimension;
%                  look at the helper visualizeLevelSet to see the options
%                  (optional, default depends on grid dimension).
%   etchShape    Shape of the etch (can be currently set to triangle or
%                rectangle)
%
%   data         Implicit surface function at t_max.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 2/9/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('../addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1.0;                  % End time.
plotSteps = 9;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 0;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 1;

%---------------------------------------------------------------------------
% Create the grid.
g.dim = 2;
g.min = -1;
g.dx = 1 / 100;
g.max = +1;
g.bdry = @addGhostExtrapolate;

g = processGrid(g);

%---------------------------------------------------------------------------
% Assign velocities
if(nargin < 1)
    etchShape =  'triangle';
end

switch(etchShape)
    case 'rectangle'
        % Velocity straight down - makes a rectangular trench
        vel{1,1} = zeros(size(g.xs{1,1}));
        vel{2,1} = zeros(size(g.xs{1,1}));

        if (g.dim > 2)
             % 3-dimensinal: velocity in the negative z-direction
             % between -0.25 < x < 0.25, and -0.25 < y < 0.25 & -0.25 < z
            vel{3,1} = zeros(size(g.xs{1,1}));
            vel{3,1}(g.xs{1,1}>-0.25 & g.xs{1,1} < 0.25 & g.xs{2,1}> -0.25 & g.xs{2,1} < 0.25 & g.xs{3,1}>-0.25) = -0.5;
        else
            % Velocity in the negative y-direction  
            % between -0.25 < x < 0.25, and 0 < y < 0.9 
            vel{2,1}(g.xs{1,1}>-0.25 & g.xs{1,1} < 0.25 & g.xs{2,1}>0 & g.xs{2,1} < 0.9) = -0.5;
        end
    case 'triangle'
        % Velocity straight down - makes a triangular trench
        vel{1,1} = zeros(size(g.xs{1,1}));
        vel{2,1} = zeros(size(g.xs{1,1}));

        if (g.dim > 2)
            % 3-dimensinal: velocity in the negative z-direction
            % between -0.25 < x < 0.25, and -0.25 < y < 0.25 & -0.25 < z
            vel{3,1} = zeros(size(g.xs{1,1}));
            vel{1,1}(g.xs{1,1}>-0.25 & g.xs{1,1} < 0 & g.xs{2,1}> -0.25 & g.xs{2,1} < 0.25) = 0.5;
            vel{1,1}(g.xs{1,1}>0 & g.xs{1,1} < 0.25 & g.xs{2,1}> -0.25 & g.xs{2,1} < 0.25) = -0.5;
            vel{3,1}(g.xs{1,1}>-0.25 & g.xs{1,1} < 0.25 & g.xs{2,1}> -0.25 & g.xs{2,1} < 0.25 & g.xs{3,1}>-0.25) = -0.5;
        else
            % Velocity in the negative y-direction  
            % between -0.25 < x < 0.25, and 0 < y < 0.9
            vel{1,1}(g.xs{1,1}>-0.6 & g.xs{1,1} < 0) = -0.1;
            vel{1,1}(g.xs{1,1}>= 0 & g.xs{1,1} < 0.6) = 0.1;
            %vel{2,1}(g.xs{1,1}>-0.6 & g.xs{1,1} < 0.6 & g.xs{2,1}>0) = -0.1;
            vel{2,1}(g.xs{2,1} > -0.8) = -0.5;
            vel{2,1}(g.xs{2,1}>0.75 & g.xs{1,1} > 0.25) = 0;
            vel{2,1}(g.xs{2,1}>0.75 & g.xs{1,1} < -0.25) = 0;
            %& g.xs{1,1} < -0.25 & g.xs{1,1} > 0.25) = 0;
            vel{1,1}(g.xs{2,1}>0.75 & g.xs{1,1} < -0.25) = 0;
            vel{1,1}(g.xs{2,1}>0.75 & g.xs{1,1} > 0.25) = 0;
        end
end

%---------------------------------------------------------------------------
% What kind of display?
if(nargin < 2)
  switch(g.dim)
   case 1
    displayType = 'plot';
   case 2
    displayType = 'contour';    
   case 3
    displayType = 'surface';
   otherwise
    error('Default display type undefined for dimension %d', g.dim);
  end
end

%---------------------------------------------------------------------------
% Create initial conditions (a rectangle).

widths = 0.8*(g.max-g.min);
center = 0.8*(g.max+g.min);
data = shapeRectangleByCenter(g,center,widths);

data0 = data;

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'low';
end

% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.velocity = @velocityFunc;
schemeData.grid = g;

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  schemeData.derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  schemeData.derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  schemeData.derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  schemeData.derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Initialize Display
f = figure;

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);

hold on;
if(g.dim > 1)
  %axis(g.axis);
  %daspect([ 1 1 1 ]);
  camlight
end

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tPlot) ];
  
  % Take a timestep.
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

  % Get correct figure, and remember its current view.
  figure(f);
  [ figure_az, figure_el ] = view;

  % Delete last visualization if necessary.
  if(deleteLastPlot)
    delete(h);
  end

  % Move to next subplot if necessary.
  if(useSubplots)
    plotNum = plotNum + 1;
    subplot(rows, cols, plotNum);
  end

  % Create new visualization.
  h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);

  % Restore view.
  view(figure_az, figure_el);
  camlight
  
end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);

function out = velocityFunc(
