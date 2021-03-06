function [ data, g, width,height ] = etchingWithOptions(g,options,current,i)

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;
%---------------------------------------------------------------------------

% Assign velocities
if(nargin < 1)
    [g, options] = createOptionsDefault();
    current = 1; i = 1;
end

% initialize mask and data
[data, mask, mask_vel, init_level] = initializeGeom(g,options);

vel = getVelocities(g,mask_vel,init_level,options, current,i);
                        
% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.velocity = vel;
schemeData.grid = g;

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');
                          
% Choose approximations at appropriate level of accuracy.
switch(options.accuracy)
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

if(options.singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Initialize Display
if(options.doDisplay)
    f = figure;
end

% Set up subplot parameters if necessary.
if(options.useSubplots)
  rows = ceil(sqrt(options.plotSteps));
  cols = ceil(options.plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

if(options.doDisplay)
    visualizeEtch(g, data, mask, options, [ 't = ' num2str(options.t0) ]);
    hold on;
end

if(g.dim > 2)
  %axis(g.axis);
  %daspect([ 1 1 1 ]);
  camlight
end

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = options.t0;
startTime = cputime;
while(options.tMax - tNow > small * options.tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(options.tMax, tNow + options.tPlot) ];
  
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
  
  
  if(options.pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

%  Get correct figure, and remember its current view.
if (options.doDisplay)
  figure(f);
  [ figure_az, figure_el ] = view;

  % Delete last visualization if necessary.
  if(options.deleteLastPlot)
    delete(h);
  end

  % Move to next subplot if necessary.
  if(options.useSubplots)
    plotNum = plotNum + 1;
    subplot(rows, cols, plotNum);
  end

  % Create new visualization.
  if(options.doDisplay)
    visualizeEtch(g, data, mask, options, [ 't = ' num2str(tNow) ]);
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

top_ind = top_ind(end);

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

function vel = getVelocities(g,mask,init_level,options,current,i)

switch(options.etchShape)
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
if(options.doMask)
 vel{1,1}(mask >= 0) = 0; vel{2,1}(mask >= 0) = 0;
end

if(g.dim > 2) 
    vel{3,1}(mask > 0) = 0;  
end

end


