function [ data, g, mask, curv] = changeBfwdModel(options,data0In,gIn,maskIn)

%---------------------------------------------------------------------------
% function [ data, g, mask, curv] = changeBfwdModel
% function [ data, g, mask, curv] = changeBfwdModel(options)
% function [ data, g, mask, curv] = changeBfwdModel(options,data0In,gIn,maskIn)
%
% Forward model
% Step1 - Compressible model (see compressModel.m)
% Step2 - Increment curvature/pressure and run constant curvature model
%         (see constCurvModel.m) as long as the level set is moving forward
%
% Need to set options.dc  - curvature increment.
% If input data is supplied, options.a (pressure corresponding to input
% data) should be set as well.
%
% Parameters:
%   'options' - structure set with setOptions() function; if empty it will be set 
%             to defaults.
%   data0In, gIn, maskIn - one can optionally provide initial conditions that are 
%            final result of the previous run. 
%  
% Output:
%   data         Implicit surface function at time options.tMax.
%   g            Grid structure on which data was computed.
%   mask         Different geometries are represented by some of the volume
%   being masked
%   curv         curvature data, size of array same as data
%---------------------------------------------------------------------------

% Make sure we can see the kernel m-files.
run('../addPathToKernel');

if( (nargin == 0) ||  ~isstruct(options ))
    disp('Options set to default values');
    options = setOptions;
end;

[AvailGeom, MaskGeom] = setAvailGeom;

fid = fopen(options.OutFileName,'w');
%get the outfile directory, data will be saved in this directory
pathstr = fileparts(options.OutFileName);
if(isempty(pathstr)) 
    pathstr = '.'; %correction for current directory
end

startTime = cputime;

if( options.doDisplay )
    fprintf(fid,'Time Start %s\n',datestr(now));
    fprintf(fid,'Forward Model Level Set Method simulation\n');
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

% get initial volume of the invading fluid 
%volume0 = size(find(data0 < 0 ));
%volume0 = volume0(1);

% get max. possible volume of the invading fluid
volume_max = size(find(mask < 0 ));
volume_max = volume_max(1)

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
elseif(doPlot)
  warning('Cannot create final plot in dimensions other than 2D');
end

step = 1;
if( nargin < 2)
    % Input data not provided, do step1
    % Figure out entry condition, i.e. entry pressure value
    if( (options.ICplane == 'x') && (g.dim == 2))
        tmp = size(find(data0(1,:) < 0));
        tmp = tmp(1) *tmp(2); %initial slab of voxels on entry
    elseif( (options.ICplane == 'y') && (g.dim == 2))
        tmp = size(find(data0(:,1) < 0));
        tmp = tmp(1) *tmp(2); %initial slab of voxels on entry    
    else disp('Fix entry condition.');
    end

    Rt = 0.5 * tmp;
    %volume_entry = 0.5 * (Rt*Rt*pi); % area of entry half-circle
    %options.Vmax_frac = volume_entry / volume_max; %corresponding volume fraction
    %fprintf(fid,'\noptions.Vmax_frac reset to %g\n',options.Vmax_frac);
    
    b0 = options.b;
    c0 = ((g.dim - 1)/(Rt*g.dx(1))); %breakthrough curvature
    a0 = b0*c0; %breakthrough pressure
    fprintf(fid,'\nEntry Rt %g, entry curvature ((g.dim-1)/Rt) = %g',Rt*g.dx(1),c0);
    
    options1 = setOptions(options,'b',b0,'a',a0,'doDisplay',0,'doSave',0);
    fprintf(fid,'\nStep %d compressible model',step);
    fprintf(fid,'\nb0 = %g (kept the same), a0 = %g (calc from entry curv.) Vmax_frac %g\n',b0,a0,options1.Vmax_frac);
    
    %Run compressible model to initialize the data.
    fclose(fid);
    [data, g, mask, curv]= compressModel(options1,data,g,mask);
    fid = fopen(options.OutFileName,'a');
    
    [k_avg(step) k_min(step) k_max(step)] = displayCurvature(data,g,mask,curv,fid,0);
    
     %output current volume of the fluid occupied part 
     volume = size(find(data < 0 ));
     volume = volume(1)
     fprintf(fid,'Invading fluid volume\n\tV(tMax) = %g, fraction total %g\n',volume,volume/volume_max);
     fprintf(fid,'Maximal void space volume\n\tVolMax %g\n',volume_max);
     actual_vmax = options.Vmax_frac*volume_max;
     fprintf(fid,'Maximal volume in a(t) model\n\tVmax = %g\n',actual_vmax);
 
     %output current a
     rel_change = (volume - actual_vmax)/actual_vmax;
     norm_speed =  a0 * exp(-options.f * rel_change);
     fprintf(fid,'Normal speed\n\ta(tMax) = %g\n',norm_speed);
     fprintf(fid,'Curvature term\n\tb*k_avg = %g\n',options.b*k_avg(step));    
 
     a0 = norm_speed; options = setOptions(options,'a',a0);
     fprintf(fid,'Pressure reset to normal_speed.\n');
     b(step) = b0; 
     c0 = a0/b0;
     curvconst(step) = c0;
    
     if( doPlot && g.dim== 2 )
        if ( ceil(step/2) == step/2 )
            [ garbage, hF ] = contour(g.xs{1}, g.xs{2}, data, lev, 'r-'); hold on;
        else
            [ garbage, hF ] = contour(g.xs{1}, g.xs{2}, data, lev, 'g-'); hold on;
        end

        set(hF, 'LineWidth', 2.0);
        axis equal
        axis(g.axis);
     end
    
    if(options.doSave)
        fname = sprintf('%s/data_step%d',pathstr,step);
        save(fname,'data');
    end
else
    fprintf(fid,'Step 1 skipped, data initialized from input data.\n');
    b0 = options.b;
    b(step) = b0;
    a0 = options.a;
    c0 = a0/b0;
    curvconst(step) = c0;
    curv = curvatureSecond(g,data);
    [k_avg(step) k_min(step) k_max(step)] = displayCurvature(data,g,mask,curv,fid,0);
end

forward_flag = 1;

volume_old = size(find(data < 0 )); volume_old = volume_old(1);

bmin =b0/3.0; %arbitrary as well; it's set so the run is limited
while( (b0 >= bmin) && forward_flag )
    step = step + 1;
    c0 = c0 + options.dc;
    curvconst(step) = c0;
    b0 = a0/c0;
    b(step) = b0;
    fprintf(fid,'\nStep %d const.curv. model \nb %g a %g\n',step,b0,a0);
    options1 = setOptions(options,'b',b0,'doDisplay',0,'doSave',0);
    fclose(fid);
    [data, g, mask, curv] = constCurvModel(options1,data,g,mask);
    fid = fopen(options.OutFileName,'a'); 
    
    volume = size(find(data < 0 )); volume = volume(1);
    % Process/display results 
    if(options.doDisplay)
         %compute current curvature
         [k_avg(step) k_min(step) k_max(step)] = displayCurvature(data,g,mask,curv,fid,0);

         %output current volume of the fluid occupied part 
         fprintf(fid,'Invading fluid volume\n\tV(tMax) = %g, fraction total %g\n',volume,volume/volume_max);     
    end
    
    if( volume < volume_old || volume == 0 )
        forward_flag = 0;
        fprintf(fid,'\nLevel set started moving backwards.');
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
    end
    
    % I think this hit_end business works really only for throat/tube
      % geometry, for Sphere2d is strange for sure
    if( (options.ICplane == 'x' ) && (g.dim == 2) )
      hit_end = size( find(data(g.N(1),:) < 0) ); 
    elseif( (options.ICplane == 'y' ) && (g.dim == 2) )
      hit_end = size( find(data(:,g.N(2)) < 0) );
    else
      hit_end = [0 0];
      disp('Fix hit_end condition');
    end

    hit_end = hit_end(1)*hit_end(2);
    if (hit_end >= 2 )  %level set hit the opposite boundary
          fprintf(fid,'\nLevel set hit the open boundary.');
          bmin = bmin - 1; % this effectively stops the iteration
    end
end 
 
%plot curvature vs. expected curvature
plotKavg(k_max,k_avg,curvconst);

if(options.doSave)
    fname = sprintf('%s/k_avg',pathstr);
    save(fname,'k_avg');
    fname = sprintf('%s/k_max',pathstr);
    save(fname,'k_max');
    fname = sprintf('%s/curvconst',pathstr);
    save(fname,'curvconst');
    fname = sprintf('%s/b',pathstr);
    save(fname,'b');
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
