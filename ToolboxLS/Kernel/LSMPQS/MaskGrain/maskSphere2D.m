function  [mask,g] = maskSphere2D(g)

% function  [mask,g] = maskSphere2D_new(g)
% Create mask that corresponds to void space btw 3 equal spheres.
% Spheres of radius 1.0 aee created (so that calculations are
% automatically "normalized" by sphere radius)

n = 3;  % number of masking spheres

maskCenter(:,1) = [0; 1.7; 0.0; 0.0 ];
maskRadius(1) = 1.0;

maskCenter(:,2) = [0; -0.6; 0.0; 0.0 ];
maskRadius(2) = 1.0;

maskCenter(:,3) = [1.75; 0.5; 0.0; 0.0 ];
maskRadius(3) = 1.0;

% The moving set can be anywhere outside the masked region.
%  '-' sign is to get the complement of the masked region.
mask = - shapeSphere(g, maskCenter(:,1), maskRadius(1));
for i=2:n
   mask = max(mask,-shapeSphere(g, maskCenter(:,i), maskRadius(i)));
end