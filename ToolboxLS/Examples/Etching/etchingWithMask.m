function [ data, g, width,height ] = etchingWithMask(accuracy, displayType, etchShape,current,i)
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
tMax = 2;                  % End time.
plotSteps = 9;               % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% To plot or not to plot
doDisplay = 1;

% Pause after each plot?
pauseAfterPlot = 0;

% Delete previous plot before showing next?
deleteLastPlot = 0;

% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 1;

%---------------------------------------------------------------------------
% Create the grid.
g.dim = 2;
g.min = [-1; -2];
g.dx = 1 / 100;
g.max = [1;2];
g.bdry = @addGhostExtrapolate;

g = processGrid(g);
                          
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

widths = (g.max-g.min)+2*g.dx;
center = g.min + 0.5*(g.max - g.min);
center(end) = center(end) - 0.2;
data = shapeRectangleByCenter(g,center,widths);

%---------------------------------------------------------------------------
% Create mask 
%remove Mask for block copolymers by setting half_width and thickness equal
%to zero
half_width =0.1;
thickness = 6*g.dx;
init_level = center(end)+widths(end)/2;

mask_widths = [(g.max(1:end-1) - g.min(1:end-1) + 2*g.dx(1)); thickness];
mask_centers = [(g.min(1:end-1) + 0.5*(g.max(1:end-1) - g.min(1:end-1))); init_level];

mask = shapeRectangleByCenter(g,mask_centers,mask_widths);
mask_vel = shapeRectangleByCenter(g,mask_centers,[mask_widths(1);mask_widths(2) + 2*g.dx]);

hole_centers = [(g.min(1:end-1) + 0.5*(g.max(1:end-1) - g.min(1:end-1))); init_level];
hole_widths = [half_width*2*(1:g.dim-1); thickness+2*g.dx];
hole = shapeRectangleByCenter(g,hole_centers,hole_widths);
hole_vel = shapeRectangleByCenter(g,hole_centers,[hole_widths(1); hole_widths(2) + 2*g.dx]);

mask = shapeDifference(mask,hole);
mask_vel = shapeDifference(mask_vel, hole_vel);

% The moving set can be anywhere outside the masked region.
%   So what we really need is the complement of the masked region.
mask = -mask;
mask_vel = -mask_vel;

% Need to ensure that the initial conditions satisfy the mask.
data_new = max(data, mask);

data = data_new;
%---------------------------------------------------------------------------
% Create layers
colors = ['b'];
[startMatInd, endMatInd, centMatInd, map] = getMaterialMap(g);
for i = 1:1
    layer_boundaries(i,:) = [endMatInd(1,i),startMatInd(1,i),];
end
layer_boundaries = flipud(layer_boundaries)
 assignin('base', 'layerboundaries', layer_boundaries)
%---------------------------------------------------------------------------

% Assign velocities
if(nargin < 1)
    etchShape =  'rectangle';
    current = 1; i = 1;
end

vel = getVelocities(g,mask_vel,init_level,etchShape, current,i);
                        
if(nargin < 1)
  accuracy = 'high';
end

% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.velocity = vel;
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
% Initialize Display%
if(doDisplay)
f = figure;

%Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

if(doDisplay)
    visualizeEtch(g, data, mask, layer_boundaries, colors, [ 't = ' num2str(t0) ]);
end

hold on;
if(g.dim > 2)
  %axis(g.axis);
  %daspect([ 1 1 1 ]);
  camlight
end
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
  [ t , y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

  %if (options.doReinit)
     %reinitialization - make data be a signed distance function
     disp('Reinitializing');
     data = signedDistanceIterative(g,data);
  %end
  
  
  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

%  Get correct figure, and remember its current view.
if (doDisplay)
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
  if(doDisplay)
    visualizeEtch(g, data, mask, layer_boundaries, colors, [ 't = ' num2str(tNow) ]);
  end
  % Restore view.
  view(figure_az, figure_el);
  %camlight
end 
end

[width, height] = getEtchDims(g, data, mask, init_level);
endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);
end

function [bottomWidth, centerHeight] = getEtchDims(g, data, mask, init_level)

intfcPoints = isNearInterface(data);
w = zeros(size(data,1),1);
h = zeros(size(data,2),1);
intfcPoints(mask >= 0) = 0;

% locate index of initial level
top_ind = find(abs(g.xs{2,1}(1,:)-init_level) < g.dx(1));

center = round(size(data,1)/2);
for i = 1:size(data,1)
    top_end = 0; bottom_end = 0;
    bottom_ind = 0;
    ind = find(intfcPoints(i,:));
    
    if(size(ind,2) == 2)
        bottom_ind = ind(2);
        bottom_end = data(i,ind(2));
        top_end = 0;
    end
    h(i) = (top_ind  - bottom_ind)*g.dx(1) + top_end - bottom_end;

end
centerHeight = h(center);

bot_trench = 1;

for j = 1:size(data,2)
    ind = find(intfcPoints(:,j));
    if(size(ind,1)>0)
        bot_trench = j;
        break;
    end
end

for j = 1:size(data,2)
    left_end = 0; right_end = 0;
    left_ind = 0; right_ind = 0;
    ind = find(intfcPoints(:,j));
    
    if(size(ind,1) == 4)
        left_end = data(ind(2),j); right_end = data(ind(3),j);
        left_ind = ind(2); right_ind = ind(3);
    end
    
    if(size(ind,1) > 4)
        left_end = data(ind(2),j); right_end = data(ind(end-1),j);
        left_ind = ind(2); right_ind = ind(end-1);
        %bot_trench = j;
    end
    w(j) = (right_ind - left_ind)* g.dx(1)+right_end - left_end;
end

bottomWidth = w(bot_trench);

end

function vel = getVelocities(g,mask,init_level,etchShape,current,i)

switch(etchShape)
    case 'rectangle'
        % Velocity straight down - makes a rectangular trench
        vel{1,1} = zeros(size(g.xs{1,1}));
        vel{2,1} = zeros(size(g.xs{1,1}));

        if (g.dim > 2)
             % 3-dimensinal: velocity in the negative z-direction
            vel{3,1} = zeros(size(g.xs{1,1}));
            vel{3,1}(:,:,:) = -0.5;
        else
            % Velocity in the negative y-direction  
             vel{2,1}(:,:) = -0.5;
        end
    case 'triangle'
        % Velocity straight down - makes a triangular trench
        vel{1,1} = zeros(size(g.xs{1,1}));
        vel{2,1} = zeros(size(g.xs{1,1}));

        if (g.dim > 2)
            % 3-dimensional: velocity in the negative z-direction
            vel{3,1} = zeros(size(g.xs{1,1}));
            vel{1,1}(g.xs{1,1} < 0) = 0.5;
            vel{1,1}(g.xs{1,1} > 0) = -0.5;
            vel{3,1}(:,:,:) = -0.5;
        else
            % Velocity in the negative y-direction  
            vel{1,1}(g.xs{1,1} < 0) = -0.1;
            vel{1,1}(g.xs{1,1} > 0) = 0.1;
            vel{2,1}(:,:) = -0.5;
        end
    case 'varyAlpha'
            % Velocity straight down - makes a triangular trench
        vel{1,1} = zeros(size(g.xs{1,1}));
        vel{2,1} = zeros(size(g.xs{1,1}));
        alpha = 0.5;
        vel_d = -0.5;
        if (g.dim > 2)
            % 3-dimensinal: velocity in the negative z-direction
            vel{3,1} = zeros(size(g.xs{1,1}));
            vel{1,1}(g.xs{1,1} < 0) = -alpha*vel_d;
            vel{1,1}(g.xs{1,1} > 0) = alpha*vel_d;
            vel{3,1}(:,:,:) = vel_d;                          
        else
            % 2-dimensional: Velocity in the negative y-direction  
            % between -0.25 < x < 0.25, and 0 < y < 0.9
%              vel{1,1}(g.xs{1,1} < 0) = .9;%alpha*vel_d;
%              vel{1,1}(g.xs{1,1} > 0) = -.9;%-alpha*vel_d;
%              vel{2,1}(:,:) = -.4; vel_d;
          [vel{1,1},vel{2,1}] = computeVelocities(g,init_level,current,i);
        end
    case 'expAlpha'
            % Velocity straight down - makes a triangular trench
        vel{1,1} = zeros(size(g.xs{1,1}));
        vel{2,1} = zeros(size(g.xs{1,1}));
        alpha = 0.5;
        vel_d = -1.2;
        lambda = 0.5;
        if (g.dim > 2)
            % 3-dimensinal: velocity in the negative z-direction
            d = g.xs{3,1} - init_level;
            v_s = alpha*exp(d/lambda)*vel_d;
            vel{3,1} = zeros(size(g.xs{1,1}));
            vel{1,1}(g.xs{1,1} < 0) = v_s(g.xs{1,1} < 0);
            vel{1,1}(g.xs{1,1} > 0) = -v_s(g.xs{1,1} > 0);
            vel{3,1}(:,:,:) = vel_d;                            
        else
            % 2-dimensional: Velocity in the negative y-direction  
            d = (g.xs{2,1} - init_level);%/(init_level + 0.8);
            v_s = alpha*exp(d/lambda)*vel_d;
            vel{1,1}(g.xs{1,1} < 0) = v_s(g.xs{1,1} < 0);
            vel{1,1}(g.xs{1,1} > 0) = -v_s(g.xs{1,1} > 0);
            vel{2,1}(:,:) = vel_d;
        end   
      
end

% % Apply mask
 vel{1,1}(mask >= 0) = 0; vel{2,1}(mask >= 0) = 0;

if(g.dim > 2) 
    vel{3,1}(mask > 0) = 0;  
end

end
