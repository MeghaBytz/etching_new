function [angle median_angle mean_angle] = explore_contact_angle(infl_name1,infl_name2,do_vis)

% [angle median_angle mean_angle] = explore_contact_angle(infl_name1,infl_name2,do_vis)
% explores contact angle btw fluid1-fluid2 interface and fluid1-grain
% interface
% infl_name1 - filename for file containing the signed distance function of
%               fluid blob of interest (produced by 3DMA-Rock)
% infl_name2 - filename with signed distance function for the complimentary
%               fluid phase in the blob subvolume (produced by 3DMA-Rock)
% do_vis - set to 1 if you want to visualaize result
% angle  = contact angle calculated at vertices of contact btw
%          fluid1-fluid2 isosurface and fluid1-grain isosurface
% median_angle,mean_angle - median and mean of the array angle


if (nargin < 3) do_vis = 0; end
    
run('~/ToolboxLS-1.0/Examples/addPathToKernel.m');
addpath ~/3dma_rock/src/matlab  %for read_l2_brnfl

% data1 < 0 defines region1 occupied by fluid1 and similarly data2
[data1 nx ny nz] = read_l2_brnfl(infl_name1);
[data2 nx ny nz] = read_l2_brnfl(infl_name2);

% color will be used if do_vis = 1
if( findstr('H2O',infl_name1) )     color = [0 0.8 0]; %green
elseif( findstr('OIL',infl_name1) ) color = [0.8 0 0]; %red
else                                color = 'y';       %yellow
end

isolevel = 0; % indicates level set of the interface ni data1 and data2
[angle median_angle mean_angle] = explore_contact_angle_main(data1,data2,nx,ny,nz,isolevel,do_vis,color);
