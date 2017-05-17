function  [mask,g] = maskThroat2Dsmall(g)

% function  [mask,g] = maskThroat2D(g)
% Create mask that corresponds to a throat in between 2 equal spheres.
% Spheres have radius 1.0 (so that calculations are
% automatically "normalized" by sphere radius)
% Gap btw the spheres (esentially throat size) can be changed below.
% Grid 'g' should be centered at 0, something like [-0.5,0.5]

n = 2;  % number of masking spheres
R = 1.0;
gap = 0.2;

maskCenter(:,1) = [0; -gap/2 - R; 0.0; 0.0 ];
maskRadius(1) = R;

maskCenter(:,2) = [0; gap/2 + R; 0.0; 0.0 ];
maskRadius(2) = R;

% The moving set can be anywhere outside the masked region.
%  '-' sign is to get the complement of the masked region.
mask = - shapeSphere(g, maskCenter(:,1), maskRadius(1));
mask = max(mask,-shapeSphere(g, maskCenter(:,2), maskRadius(2)));
    