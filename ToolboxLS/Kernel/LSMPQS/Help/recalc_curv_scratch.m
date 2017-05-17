function [k_max k_avg] = recalc_curv_scratch(pathstr,num)

% function  [k_max k_avg] = recalc_curv_scratch(pathstr,num)

addpath('/ices/masha/ToolboxLS-1.0/Kernel/Masha/Help');

fname = sprintf('%s/mask',pathstr);
mask = read_double_array(fname);
fname = sprintf('%s/grid',pathstr);
g = read_grid_file(fname);

fname = sprintf('%s/curvconst',pathstr);
curvconst = read_double_array1d(fname);

% form Toolbox style grid structure from provided data
g1.min = g.x_lo(1:2);
g1.max = g.x_hi(1:2);
g1.dim = 2;
%g1.dx = g.dx(1:2);
g1.N = g.n(1:2);

g1 = processGrid(g1);

fid = fopen('garbage','w');

for(i=1:num)
    fname = sprintf('%s/data_step%d',pathstr,i);
    data = read_double_array(fname);

    [curv grad_norm_data] = curvatureSecond(g1,data);
    [k_avg(i) k_min(i) k_max(i)] = displayCurvature(data,g1,mask,curv,fid,0);
end

fclose(fid);
system('/bin/rm garbage');

plotKavg(k_max,k_avg,curvconst);