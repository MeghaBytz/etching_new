function  [mask,g] = maskThroatBump1(g)

% function  [mask,g] = maskThroatBump1(g)

n = 2;  % number of masking spheres
R = 1.0;
gap = 0.3;

maskCenter(:,1) = [0; -gap/2 - R; 0.0; 0.0 ];
maskRadius(1) = R;

maskCenter(:,2) = [0; gap/2 + R; 0.0; 0.0 ];
maskRadius(2) = R;

% The moving set can be anywhere outside the masked region.
%  '-' sign is to get the complement of the masked region.
mask = - shapeSphere(g, maskCenter(:,1), maskRadius(1));
mask = max(mask,-shapeSphere(g, maskCenter(:,2), maskRadius(2)));

maskCenter(:,3) = [0.3; -0.25; 0.0; 0.0 ];
maskRadius(3) = 0.1;

mask = max(mask,-shapeSphere(g, maskCenter(:,3), maskRadius(3)));

    