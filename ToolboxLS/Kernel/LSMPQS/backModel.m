function [ data, g, mask, curv] = backModel(options,data0In,gIn,maskIn,amin)

%--------------------------------------------------------------------------
% function [ data, g, mask, curv] = backModel(options,data0In,gIn,maskIn)
%
% Level Set Method Progressive Quasistatic Algorithm for Imbibition.
% Step1 - read in the data (end point of drainage).
%         options.a and options.b should correspond to teh drainage end point
% At each further step - Decrement curvature/pressure and run constant 
%         curvature model (see constCurvModel.m) as long as the level set 
%         is moving backward (region data < 0 is getting smaller).
%
% Pressure increment is options.dc should be negative.
% Set options.a (pressure) to correspond to input data.
%
% Parameters:
%   options - structure set with setOptions() function; if empty it will be set 
%             to defaults.
%   data0In, gIn, maskIn  - initial conditions that are 
%            final result of a previous drainage run. 
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


if( (nargin == 0) ||  ~isstruct( options ))
    disp('Options set to default values');
    options = setOptions;
end;

[AvailGeom, MaskGeom] = setAvailGeom;

fid = fopen(options.OutFileName,'w');
%get the outfile directory, data will be saved in this directory
pathstr = fileparts(options.OutFileName);
if(isempty(pathstr)) 
    pathstr = '.'; %correction for the current directory
end

startTime = cputime;

if( options.doDisplay )
    fprintf(fid,'Time Start %s\n',datestr(now));
    fprintf(fid,'Progressive Quasi-static Model Level Set Method simulation\n');
    printOptions(options,fid);
end

%---------------------------------------------------------------------------

if( options.doDisplay )
    % Note that there are some additional display options below that can be
    % changed only within the file.
    doPlot = 1;                  % Produce plot after each step?
    %useSubplots = 0;
else
    doPlot = 0;
end

if( options.doDisplay )
    fprintf(fid,'Simulation continues until given tMax is reached\n'); 
    fprintf(fid,'or max.abs.error for data(:,tFinal)-data(:,tFinal-tPlot');
    fprintf(fid,'is less than epsStop.\n');
    fprintf(fid,'Simulation will also stop if max.abs.error starts increasing.\n');
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
    error('Input data not provided.');
else
    data0 = data0In;
    data  = data0In;
    g     = gIn;
    mask  = maskIn;
    fprintf(fid,'Data initialized from user provided input.\n');
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

% get initial volume of the Retreating fluid 
%volume0 = size(find(data0 < 0 ));
%volume0 = volume0(1);

% get max. possible volume of the Retreating fluid
volume_max = size(find(mask < 0 ));
volume_max = volume_max(1);

if(doPlot && g.dim == 2 )
  % Display initial set, mask
  figure;
  lev = [ level level ];
  [ garbage, hI ] = contour(g.xs{1}, g.xs{2}, data0, lev,'b-');
  hold on;
  [ garbage, hM ] = contourf(g.xs{1}, g.xs{2}, mask, lev, 'k-');colormap gray;
  
  hs = [ hI; hM ];

  set(hs, 'LineWidth', 2.0);

  legend([ hI(1), hM(1) ], ...
         {'initial','mask'}, 2)
  axis equal
  axis(g.axis);
  pause(1);
elseif(doPlot)
   disp('Cannot create final plot in dimensions other than 2D');
end

if( options.dc >= 0) 
    fprintf('\noptions.dc=%g is positive. Changed to %g.\n',options.dc,...
                                                           -options.dc);
    options.dc = -options.dc;
end

da = options.dc*options.b;
fprintf(fid,'\nPressure decrement da = %d\n',da);

if(options.doSave)
    fname = sprintf('%s/data_init',pathstr);
    save(fname,'data0');
    fname = sprintf('%s/mask',pathstr);
    save(fname,'mask');
    fname = sprintf('%s/grid',pathstr);
    save(fname,'g');
end

step = 1; %input data
a0 = options.a;
a(step) = a0;
curvconst(step) = options.a/options.b;
curv = curvatureSecond(g,data);
[k_avg(step) k_min(step) k_max(step)] = displayCurvature(data,g,mask,curv,fid,0);

backward_flag = 1;

%data < 0 is where the non-wetting fluid is
volume_old = size(find(data < 0 )); volume_old = volume_old(1);

if (nargin < 5 )
    amin = 0.005; %arbitrary; it's set so the run is limited
end

while( (a0 >= amin) && backward_flag )
    step = step + 1;
    a0 = a0 + da;
    curvconst(step) = a0/options.b;
    a(step) = a0;
    fprintf(fid,'\nStep %d const.curv. model \na %g\n',step,a0);
    options1 = setOptions(options,'a',a0,'doDisplay',0,'doSave',0);
    fclose(fid);
    [data, g, mask, curv] = constCurvModel(options1,data,g,mask);
    fid = fopen(options.OutFileName,'a'); 
    
    volume = size(find(data < 0 )); volume = volume(1);
    % Process/display results 
    if(options.doDisplay)
         %compute current curvature
         [k_avg(step) k_min(step) k_max(step)] = displayCurvature(data,g,mask,curv,fid,0);

         %output current volume of the fluid occupied part 
         fprintf(fid,'Retreating fluid volume\n\tV(tMax) = %g, fraction total %g\n',volume,volume/volume_max);     
    end
    
    if( (volume > volume_old + 2) || volume == 0 )
        backward_flag = 0;
        fprintf(fid,'\nLevel set started moving forward.');
    else
        volume_old = volume;
        
        if(options.doSave)
            fname = sprintf('%s/data_step%d',pathstr,step);
            save(fname,'data');
        end
    end
    
    if( doPlot && g.dim== 2 )
        if ( ceil(step/2) == step/2 )
            [ garbage, hF ] = contour(g.xs{1}, g.xs{2}, data, lev, 'r-'); hold on;
        else
            [ garbage, hF ] = contour(g.xs{1}, g.xs{2}, data, lev, 'g-'); hold on;
        end

        set(hF, 'LineWidth', 2.0);
        axis equal
        axis(g.axis);
        pause(1);
    end
    
end 
 
%plot curvature vs. expected curvature
%plotKavg(k_max,k_avg,curvconst);

if(options.doSave)
    fname = sprintf('%s/k_avg',pathstr);
    save(fname,'k_avg');
    fname = sprintf('%s/k_max',pathstr);
    save(fname,'k_max');
    fname = sprintf('%s/curvconst',pathstr);
    save(fname,'curvconst');
    fname = sprintf('%s/a',pathstr);
    save(fname,'a');
end

%---------------------------------------------------------------------------
% Time and memory summary
endTime = cputime;
fprintf('\nTotal execution time %g seconds', endTime - startTime);
% This is printed to file even if options.doDisplay==0
fprintf(fid,'\nTotal execution time %g seconds', endTime - startTime);

disp('Memory usage after main loop'); whos
s = whos;
[n m]= size(s); mem=0;
for i=1:n, mem=mem+s(i).bytes; end
fprintf('\nTotal memory used (after main loop) %g bytes',mem);
if (options.doDisplay)
    fprintf(fid,'\nTotal memory used (after main loop) %g bytes\n',mem);
end

%---------------------------------------------------------------------------

fclose(fid);
