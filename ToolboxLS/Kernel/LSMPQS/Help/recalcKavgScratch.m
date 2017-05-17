function [k_max k_avg] = recalcKavgScratch(pathstr,num)

% function  [k_max k_avg] = recalcKavgScratch(pathstr,num)

fname = sprintf('%s/mask.mat',pathstr);
load(fname);
fname = sprintf('%s/grid.mat',pathstr);
load(fname);

fname = sprintf('%s/curvconst.mat',pathstr);
load(fname);

mask = signedDistanceIterative(g,mask,'veryHigh',5*g.dx(1),1e-3);

fid = fopen('garbage','w');

for(i=1:num)
    fname = sprintf('%s/data_step%d.mat',pathstr,i);
    load(fname);

    [curv grad_norm_data] = curvatureSecond(g,data);
    [k_avg(i) k_min(i) k_max(i)] = displayCurvature(data,g,mask,curv,fid,0);
end

fclose(fid);
system('/bin/rm garbage');

plotKavg(k_max,k_avg,curvconst);
