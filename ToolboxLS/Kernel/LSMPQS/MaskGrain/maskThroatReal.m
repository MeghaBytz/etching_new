function [mask,g] = maskThroatReal(g)
% [mask,g] = maskThroatReal(g)
% Create mask by reading in a signed distance function (a specific 3DMA l2 burn file 
% g - grid structure
% Grid structure is determined by dimensions store in the file

fid = fopen('RealInput/real_thr_sgn_dist','r');
e  = fread(fid,1,'char');
nx = fread(fid,1,'int');
ny = fread(fid,1,'int'); 
zs = fread(fid,1,'int');
ze = fread(fid,1,'int');

% we assume zs=ze, i.e. 2D data 
mask = fread(fid,[nx ny],'single');
fclose(fid);

% Grid dimensions determined by the dimensions from file
g.dim = 2;
g.min = [1;1]*g.dx;
g.max = [nx;ny]*g.dx;
mask = g.dx*mask;