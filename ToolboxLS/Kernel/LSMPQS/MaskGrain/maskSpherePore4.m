function  [mask,g] = maskSpherePore4(g)

% function  [mask,g] = maskSphere2D_new(g)
% Create mask that corresponds to void space btw 4 equal spheres.
% Spheres of radius 1.0 aee created (so that calculations are
% automatically "normalized" by sphere radius)

n = 4;  % number of masking spheres
R = 1.0;

h = 1.0;

maskCenter(:,1) = [-h; -h; 0.0; 0.0 ];
maskRadius(1) = R;

maskCenter(:,2) = [0.96*h; -0.98*h; 0.0; 0.0 ];
maskRadius(2) = R;

maskCenter(:,3) = [-h;h; 0.0; 0.0 ];
maskRadius(3) = R;

maskCenter(:,4) = [0.96*h; 0.98*h; 0.0; 0.0 ];
maskRadius(4) = R;

% The moving set can be anywhere outside the masked region.
%  '-' sign is to get the complement of the masked region.
mask = - shapeSphere(g, maskCenter(:,1), maskRadius(1));
for i=2:n
   mask = max(mask,-shapeSphere(g, maskCenter(:,i), maskRadius(i)));
end