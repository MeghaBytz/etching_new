function [ data, g, mask, curv] = compressModelCA(options,data0In,gIn,maskIn)

%---------------------------------------------------------------------------
% function [ data, g, mask, curv] = compressModel(options)
% function [ data, g, mask, curv] = compressModel(options,data0In,gIn,maskIn)
%
% Constant curvature level set simulation
% Governing PDE is
%           phi_t(x,t) + a(t)| grad phi(x,t)| = b kappa(phi) |grad phi(x,t)|
% phi - implicit function whose interior represents invading fluid
% _t  denotes temporal derivative
% grad denotes spatial gradient
% kappa denotes curvature
% a(t) is capillary pressure modeled as a(t) = exp(-f*(V(t)- Vmax)/Vmax) where V(t) is 
% volume of the invading fluid at time t, and Vmax the maximal volume the
% fluid could occupy
% f is an arbitary constant that is related to the volume of the invading
% fluid when the equation reaches the steady state

% b is a constant that models surface tension btw invading and resident fluid

% Parameters:
%   'options' - structure set with setOptions() function; if empty it will be set 
%             to defaults.
%   data0In, gIn, maskIn - one can optionally provide initial conditions that are 
%            final result of the previous run. 
%  
% Output:
%   data         Implicit interface function (final).
%   g            Grid structure on which data was computed.
%   mask         Level set function corresponding to masking geometry.
%   curv         Curvature of 'data' (pointwise)
%
% Note: 
%   - in this application data=0 and mask=0 are level sets of interest
%   - Make sure we can see the kernel m-files.
%     (run ToolboxLS function addPathToKernel.m)
%
% Author: Masa Prodanovic, UT Austin
%--------------------------------------------------------------------------


if( ~isstruct(options ))
    disp('Options set to default values');
    options = setOptions;
end;

[AvailGeom, MaskGeom] = setAvailGeom;

fid = fopen(options.OutFileName,'a');

if( options.doDisplay )
    fprintf(fid,'Time Start %s\n',datestr(now));
    fprintf(fid,'Compressible Model Level Set Method simulation\n');
    printOptions(options,fid);
end

%---------------------------------------------------------------------------
% Integration parameters.
% Start time.
t0 = 0;                      
% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

if( options.doDisplay )
    % Note that there are some additional display options below that can be
    % changed only within the file.
    doIntermedPlot = 0;          % Produce intermediate plots?
    plotSteps = 4;               % How many intermediate plots to produce?
                                 % if doIntermedPlot = 0, this can still be used to
                                 % produce some other intermediate results

    singleStep = 0;              % Plot at each timestep (overrides options.tPlot).

    plotGradient = 0;
    
    % Period at which intermediate plots should be produced.
    if(doIntermedPlot) 
        options.tPlot = (options.tMax - t0) / (plotSteps - 1);
        disp('tPlot overwritten!');
    end

    % Pause after each plot?
    pauseAfterPlot = 1;
    % Delete previous plot before showing next?
    deleteLastPlot = 0;
    % Plot in separate subplots (set deleteLastPlot = 0 in this case)?
    useSubplots = 0;
else
    doIntermedPlot = 0;
    singleStep = 0;
end



if( options.doDisplay )
    fprintf(fid,'Simulation continues until given tMax is reached\n'); 
    fprintf(fid,'or max.abs.error for data(:,tFinal)-data(:,tFinal-tPlot');
    fprintf(fid,'is less than epsStop.\n');
    fprintf(fid,'Simulation will also stop if max.abs.error starts increasing.');
    fprintf(fid,'-------------------------------------------------------------\n');
end

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;
%---------------------------------------------------------------------------
% Use periodic boundary conditions?
periodic = 0;

% Create initial conditions, grid and mask if they are not provided
if (nargin < 2)
    % Initialize grid dimensions.
    [data, data0, g, mask] = initializeDataGrid(options,MaskGeom,periodic);
else
    data0 = data0In;
    data  = data0In;
    g     = gIn;
    mask  = maskIn;
end

if( options.doDisplay )
    fprintf(fid,'Grid spacing  dx %g\n',g.dx(1));
    fprintf(fid,'Volume limits [%g,%g]x[%g,%g]\n',g.min(1),g.max(1),g.min(2),g.max(2));
end

%---------------------------------------------------------------------------
% What kind of display?
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

%---------------------------------------------------------------------------
% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
%   Same accuracy is used by both components of motion.
switch(options.accuracy)
 case 'low'
  derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', options.accuracy);
end

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

% get initial volume of the invading fluid 
%volume0 = size(find(data0 < 0 ));
%volume0 = volume0(1);

% get max. possible volume of the invading fluid
volume_max = size(find(mask < 0 ));
volume_max = volume_max(1);

%---------------------------------------------------------------------------
% Set up basic motion in the normal direction.
normalFunc = @termNormal;
normalData.grid = g;
normalData.speed = @pressureModelVmax; % handle to the function that computes normal speed a(t)
normalData.magnitude = options.a;
normalData.derivFunc = derivFunc;
vol_factor = options.Vmax_frac;
normalData.vmax = vol_factor*volume_max; % used in pressureModelVmax
normalData.factor = options.f;   % used in pressureModelVmax

%---------------------------------------------------------------------------
% Set up curvature motion.
curvatureFunc = @termCurvature;
curvatureData.grid = g;
curvatureData.curvatureFunc = @curvatureSecond;
curvatureData.b = options.b;

%---------------------------------------------------------------------------
% Combine components of motion.
if(options.b > 0)
  % If there is a nonzero curvature contribution to speed.
  schemeFunc = @termSum;
  schemeData.innerFunc = { normalFunc; curvatureFunc };
  schemeData.innerData = { normalData; curvatureData };
else
  % Otherwise ignore curvature.
  schemeFunc = normalFunc;
  schemeData = normalData;
end

%---------------------------------------------------------------------------
% Set up data required for the mask operation.
%   Mask will be compared to vector form of data array used by integrator.
normalData.mask = mask(:);
normalData.doMask = options.doMask;
curvatureData.mask = mask(:);
curvatureData.doMask = options.doMask;

schemeData.mask = mask(:);
schemeData.doMask = options.doMask;

% Also keep track of minimum of phi over time.
%   Minimum will be kept in vector form used by integrator.
normalData.min = data(:);
normalData.doMin = options.doMin;
curvatureData.min = data(:);
curvatureData.doMin = options.doMin;
 
schemeData.min = data(:);
schemeData.doMin = options.doMin;

% Let the integrator know what function to call.
integratorOptions = odeCFLset(integratorOptions, ...
                              'postTimestep', @maskAndKeepMinCA);

if( doIntermedPlot)                         
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
    %plotLevelSetInterior(data,level,mask)

    hold on;
    contourf(g.xs{1}, g.xs{2}, mask, [level level], 'k-'); colormap gray;
    
    if(g.dim > 1)
      axis(g.axis);
      axis equal;
      %daspect([ 1 1 1 ]);
    end
end

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff) or max_abs_error is
% satisfactory or max_abs_error starts increasing
tNow = t0;
startTime = cputime;
max_abs_err = 1000.0;
go_on = 1;

while( go_on )
  
  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(options.tMax, tNow + options.tPlot) ];
  
  % Take a timestep.
  % Note - there is no operator splitting options here
  [ t y schemeData ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                             integratorOptions, schemeData);
                         
  max_abs_err_old = max_abs_err;                           
  max_abs_err = max(abs(y - y0));
  
  %fprintf('Time interval [%g,%g], max_abs error %g\n', tNow,t(end),max_abs_err);
  % This is the bare minimum printed to file even if options.doDisplay==0
  fprintf(fid,'Time interval [%g,%g], max_abs error %g\n', tNow,t(end),max_abs_err);
  
  tNow = t(end);
  
  % Get back the correctly shaped data array
  data = reshape(y,g.shape);
  data_old = reshape(y0,g.shape);

  non_empty = size( find (y < 0 ) ); non_empty = non_empty(1)*non_empty(2);
  
 % go_on = (options.tMax - tNow > small * options.tMax) & non_empty &...
  %        (max_abs_err > options.epsStop) & (max_abs_err < max_abs_err_old);
  go_on = (options.tMax - tNow > small * options.tMax) & non_empty &...
          (max_abs_err > options.epsStop);
 % if(max_abs_err == max_abs_err_old)
 %     fprintf('Possibly stuck with too small tPlot=%g, increase and rerun.',options.tPlot);
 % end
  
  if( doIntermedPlot)
      if(pauseAfterPlot)
        % Wait for last plot to be digested.
        pause(1);
      end

      % Get correct figure, and remember its current view.
      figure(f);
      figureView = view;

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
      axis equal
      axis(g.axis);
      if( useSubplots )
          hold on; contourf(g.xs{1}, g.xs{2}, mask,[level level], 'k-'); colormap gray;
      end
      %plotLevelSetInterior(data,level,mask)

      % Restore view.
      view(figureView);
  end
  
  if (options.doReinit)
     %reinitialization - make data be a signed distance function
     %disp('Reinitializing');
     data = signedDistanceIterative(g,data);
  end
end

%---------------------------------------------------------------------------
% Time and memory  summary
endTime = cputime;
%fprintf('Total execution time %g seconds', endTime - startTime);
% This is printed to file even if options.doDisplay==0
fprintf(fid,'Total execution time %g seconds\n', endTime - startTime);

%disp('Memory usage after main loop'); whos
s = whos;
[n m]= size(s); mem=0;
for i=1:n, mem=mem+s(i).bytes; end
%fprintf('Total memory used (after main loop) %g bytes',mem);
if (options.doDisplay)
    fprintf(fid,'Total memory used (after main loop) %g bytes\n',mem);
end

%---------------------------------------------------------------------------

%compute current curvature
curv= curvatureSecond(g,data);

% Process and display final results.
if(options.doDisplay)
    if(g.dim == 2)
      % Display initial set, mask, final set.
      figure;
      lev = [ level level ];
      [ garbage, h(1) ] = contour(g.xs{1}, g.xs{2}, data0, lev, 'b--');
      hold on;
      [ garbage, h(2) ] = contour(g.xs{1}, g.xs{2}, data, lev, 'r-');
      [ garbage, h(3) ] = contourf(g.xs{1}, g.xs{2}, mask, lev, 'k-'); colormap gray; 
      %[ garbage, h(4) ] = contour(g.xs{1}, g.xs{2}, data_old, lev, 'k:');
      set(h, 'LineWidth', 2.0);
      legend([ h(1), h(2), h(3) ], ...
             {'initial', 'final', 'mask'}, 2);
      %legend([ hI(1), hF(1), hT(1), hM(1) ], ...
      %       {'initial', 'final', 'plotStep before final', 'mask'}, 2);
      axis equal
      axis(g.axis);

      if (plotGradient)
        plotGradient2D(data0,g,mask);
        plotGradient2D(data,g,mask);
      end 
    else          
      warning('Currently not displaying plots in dimensions other than 2D');
    end

     k_avg = displayCurvature(data,g,mask,curv,fid);

     %max absolute error
     max_abs_err = max(max( abs(data - data_old)));
     fprintf(fid,'Max abs difference in  (data(:,%g) - data(:,%g)) everywhere\n',...
                                                           tNow,tNow-options.tPlot);
     fprintf(fid,'\t%g\n',max_abs_err);

     %output current volume of the fluid occupied part 
     volume = size(find(data < 0 ));
     volume = volume(1);
     fprintf(fid,'Invading fluid volume\n\tV(tMax) = %g, fraction total %g\n',volume,volume/volume_max);
     fprintf(fid,'Maximal void space volume\n\tVolMax %g\n',volume_max);
     fprintf(fid,'Maximal volume in a(t) model\n\tVmax = %g\n',vol_factor*volume_max);
 
     %output current a
     rel_change = (volume - 0.75*volume_max)/(0.75*volume_max);
     norm_speed =  options.a * exp(-vol_factor*rel_change);
     fprintf(fid,'Normal speed\n\ta(tMax) = %g\n',norm_speed);
     %fprintf(fid,'Curvature term\n\tb*k_avg = %g\n',b*k_avg);

     fprintf(fid,'Time End %s\n',datestr(now));
end

fclose(fid);
