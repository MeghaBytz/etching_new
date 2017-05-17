function [k_avg k_min k_max] = replotCurvature(data_fname,mask_fname,grid_fname)

% function replotCurvature(data_fname,mask_fname,grid_fname)
% Open filenames that store data, mask and grid and replot curvature.


load(data_fname)
load(mask_fname)
load(grid_fname)

[curv grad_norm_data] = curvatureSecond(g,data);
[k_avg k_min k_max] = displayCurvature(data,g,mask,curv);
    
    
    
    
    


