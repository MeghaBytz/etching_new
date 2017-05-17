function [angle median_angle mean_angle] = explore_contact_angle_ubc(infl_name1,infl_name2,do_vis,nx,ny,nz)

% [angle median_angle mean_angle] = explore_contact_angle_ubc(infl_name1,infl_name2,do_vis,nx,ny,nz)
% explores contact angle btw fluid1-fluid2 interface and fluid1-grain
% interface
% infl_name1,2 - filename for file containing fluids blob of interest
%              (these are segmented data unsigned binary char arrays!)
%              the first file is treated as the main fluid of interest
% do_vis - set to 1 if you want to visualaize result
% angle  = contact angle calculated at vertices of contact btw
% fluid1-fluid2 isosurface and fluid1-grain isosurface
% median_angle and mean_angle are simple median and mean of the array angle

if (nargin < 3) do_vis = 0; end
    
run('~/ToolboxLS-1.0/Examples/addPathToKernel.m');
addpath ~/3dma_rock/src/matlab  %for read_3dma_segfl

% data1 < 0 defines region1 occupied by fluid1 and similarly data2
% arrays read in below are 'uint8'
data1 = read_ubc_segfl(infl_name1,nx,ny,nz); % 0 where fluid1, otherwise 1
data2 = read_ubc_segfl(infl_name2,nx,ny,nz); % 0 where fluid2, otherwise 1

isolevel = 0.5; % indicates level set of the interface ni data1 and data2

% color will be used if do_vis = 1
if( findstr('H2O',infl_name1) )     color = [0 0.8 0]; %green
elseif( findstr('OIL',infl_name1) ) color = [0.8 0 0]; %red
else                                color = 'y';       %yellow
end

[angle median_angle mean_angle] = explore_contact_angle_main(data1,data2,nx,ny,nz,isolevel,do_vis,color)
