function [ data, g, mask, curv] = fwdModelCA_3d(options,data0In,gIn,maskIn)

%---------------------------------------------------------------------------
% function [ data, g, mask, curv] = fwdModel
% function [ data, g, mask, curv] = fwdModel(options)
% function [ data, g, mask, curv] = fwdModel(options,data0In,gIn,maskIn)
%
% Level Set Method Progressive Quasistatic Algorithm for Drainage.
% Step1 - Compressible model (see compressModel.m)
% Step >=2 - Increment/decrement curvature/pressure and run constant 
%         curvature model (see constCurvModel.m) as long as the 
%         level set is moving forward (region data < 0 is growing)
%
% Default pressure increment is options.dc, change as needed.
% If input data is supplied, options.a (pressure corresponding to input
% data) should be set as well.
%
% Parameters:
%   options - structure set with setOptions() function; if empty it will be set 
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
%{
if( (nargin == 0) ||  ~isstruct( options ))
    disp('Options set to default values');
    options = setOptionsCA;
end;

[AvailGeom, MaskGeom] = setAvailGeom;
%}
fid = fopen(options.OutFileName,'w');
%get the outfile directory, data will be saved in this directory
pathstr = fileparts(options.OutFileName);
if(isempty(pathstr)) 
    pathstr = '.'; %correction for the current directory
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

g = gIn;
%{
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
%}
%---------------------------------------------------------------------------
%{
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
%
% get initial volume of the invading fluid 
%volume0 = size(find(data0 < 0 ));
%volume0 = volume0(1);

% get max. possible volume of the invading fluid
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
  warning('Not set to plot in dimensions other than 2D');
end
%
data0 = data0In;

if(options.doSave)
    fname = sprintf('%s/data_init',pathstr);
    save(fname,'data0');
    fname = sprintf('%s/mask',pathstr);
    save(fname,'mask');
    fname = sprintf('%s/grid',pathstr);
    save(fname,'g');
end
%}
fname = sprintf('%s/data_init',pathstr);
save(fname,'data0In');

mask = maskIn;
%{
if( g.dim == 2 )
    opposite_side = findOppSide_3d(mask,g,data0);
    opp_size = numel(opposite_side);
end
%}
step = 1;
%{
if( nargin < 2)
    % Input data not provided, do step1
    % Figure out entry condition, i.e. entry pressure value
   
    pos = options.IC_pos;
    if( (options.ICplane == 'x' || options.ICplane == 'c')  && (g.dim == 2))
        dim = 1;
    elseif( (options.ICplane == 'y') && (g.dim == 2))
        dim = 2;
    else error('Fix entry condition.');
    end

    % count voxels in connected components of the invading fluid at entry
    
    [ind comp sz_comp] = findConnComp(data0,dim,pos);
    %Minimal pressure to enter is corresponding to the biggest throat
    Rt = 0.5 * max(sz_comp);

    %MODIFIED
    a_entry = 0.05 * ((g.dim - 1)/(Rt*g.dx(1))); %breakthrough pressure
    %a_entry = 0.4; % forcing a specific value
    %END
    
    fprintf(fid,'\nEntry Rt %g, entry pressure a0 = b * ((g.dim-1)/Rt) = %g',Rt*g.dx(1),a_entry);
    
    fprintf(fid,'\nStep %d compressible model',step);
    a0 = 0.75*a_entry;  
    options1 = setOptionsCA(options,'a',a0,'doDisplay',0,'doSave',0, 'doreinit',0);
    fprintf(fid,'\na0 = %g (reduced), Vmax_frac %g\n',a0,options1.Vmax_frac);
    
    %Run compressible model to initialize the data.
    fclose(fid);
    [data, g, mask, curv]= compressModelCA(options1,data,g,mask);
    fid = fopen(options.OutFileName,'a');
    
    [k_avg(step) k_min(step) k_max(step)] = displayCurvature(data,g,mask,curv,fid,0);
    
     %output current volume of the fluid occupied part 
     volume = size(find(data < 0 ));
     volume = volume(1);
     fprintf(fid,'Invading fluid volume\n\tV(tMax) = %g, fraction total %g\n',volume,volume/volume_max);
     fprintf(fid,'Maximal void space volume\n\tVolMax %g\n',volume_max);
     actual_vmax = options.Vmax_frac*volume_max;
     fprintf(fid,'Maximal volume in a(t) model\n\tVmax = %g\n',actual_vmax);
 
     %output current a
     rel_change = (volume - actual_vmax)/actual_vmax;
     norm_speed =  a0 * exp(-options.f * rel_change);
     fprintf(fid,'Normal speed\n\ta(tMax) = %g\n',norm_speed);
     fprintf(fid,'Curvature term\n\tb*k_avg = %g\n',options.dx*k_avg(step));      %MODIFIED: options.b-options.dx
     
     curvconst(step) = norm_speed;
     a0 = max(norm_speed,a_entry); fprintf(fid,'Pressure reset to max(normal_speed,entry pressure).\n');
     a(step) = a0;
    
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
    
    if(options.doSave)
        fname = sprintf('%s/data_step%d',pathstr,step);
        save(fname,'data');
    end
else
    fprintf(fid,'Step 1 skipped, data initialized from input data.\n');
    %a0 = options.a;           
    %a(step) = a0;
    %curvconst(step) = options.a(step)/options.b(step);
    curv = curvatureSecond(g,data);
    [k_avg(step) k_min(step) k_max(step)] = displayCurvature(data,g,mask,curv,fid,0);
    
end
%}
forward_flag = 1;
data = data0In;

volume_old = size(find(data < 0 )); volume_old = volume_old(1);
a0 = 0.05;
amax = 50*a0; %arbitrary as well; it's set so the run is limited
while( (a0 <= amax) && forward_flag )
    step = step + 1;
    %da = options.dc*options.dx;
    a0 = a0 + options.da;
    %k0 = k0 + options.dc;
    [a,b,vel] = abplot_new_1(mask,options.dx,a0,options.ca, g);

    options1 = setOptionsCA(options,'a',a,'b',b,'velocity',vel,'doDisplay',0,'doSave',0);
    fclose(fid);
    [data, g, mask, curv] = constCurvModelCA(options1,data,g,mask);
    fid = fopen(options.OutFileName,'a'); 
    
    volume = size(find(data < 0 )); volume = volume(1);
    % Process/display results 
    if(options.doDisplay)
         %compute current curvature
         [k_avg(step) k_min(step) k_max(step)] = displayCurvature(data,g,mask,curv,fid,0);

         %output current volume of the fluid occupied part 
         fprintf(fid,'Invading fluid volume\n\tV(tMax) = %g, fraction total %g\n',volume,volume/volume_max);     
    end
    
    if( (volume + 5< volume_old) || volume == 0 )
        forward_flag = 1;
        fprintf(fid,'\nLevel set started moving backwards.');
    else
        volume_old = volume;
        
        if(options.doSave)
            fname = sprintf('%s/data_step%d',pathstr,step);
            save(fname,'data');
        end
    end
    %{
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
    
    opp_side_invaded = numel( find( data(opposite_side) <= 0 ) );
    if( ( options.stopTouch && opp_side_invaded ) ||...
        ( opp_side_invaded == opp_size ) )
        disp('Level set hit/invaded the open boundary opposite of the IC.');
        a0 = amax + 1; % this effectively stops iteration
    end
    %}
    
    end_slice(:,:) = mask(g.N(1),:,:);
    [mask_x mask_y] = find(end_slice<0);
    
    data_slice(:,:) = data(g.N(1),:,:);
    [data_x data_y] = find(data_slice<0);
    
    mask_end = [mask_x,mask_y];
    data_end = [data_x,data_y];
    
    if(numel(intersect(mask_end,data_end)))
        a0 = amax + 1;
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
%
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
