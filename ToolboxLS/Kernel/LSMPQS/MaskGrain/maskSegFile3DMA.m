function [mask,g] = maskSegFile3DMA(g,fname)

% [mask,g] = maskSegFile3DMA(g,fname)
% Create mask by reading in a 3DMA segmented file
% g - grid structure (needs to have g.dx set to desired voxel length)
%     Grid structure dimensions will be set from the file
% fname = segmented file name

[mask nx ny nz] = readSegfl3DMA(fname);

% note that segmented data values are 0 in pore space and we need them
% to be negative 
mask( find (mask == 0 ) ) = -1;
mask = g.dx*mask; % some scaling

if( nz == 1)
    % i want things fliped for visualization
    % this is I guess just my preference to flip things 
    % because Matlab shows things in a specific way
    % size of array is [ny nx] now
    mask = mask';
end

% Grid dimensions determined by the file dimensions
if( nz == 1)
    g.dim = 2;
    g.min = [1;1]*g.dx;
    g.max = [ny;nx]*g.dx;
else
    g.dim = 3;
    g.min = [1;1;1]*g.dx;
    g.max = [nx;ny;nz]*g.dx;
end



