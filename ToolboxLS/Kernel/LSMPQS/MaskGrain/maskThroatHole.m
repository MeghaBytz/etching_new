function  [mask,g] = maskThroatHole(g)

% function  [mask,g] = maskThroatHole(g)

n = 2;  % number of masking spheres
R = 1.0;
gap = 0.3;

maskCenter(:,1) = [0; -gap/2 - R; 0.0; 0.0 ];
maskRadius(1) = R;

maskCenter(:,2) = [0; gap/2 + R; 0.0; 0.0 ];
maskRadius(2) = R;

maskCenter(:,3) = [-0.3; -0.15; 0.0; 0.0 ];
maskRadius(3) = 0.1;

mask1 = shapeSphere(g, maskCenter(:,3), maskRadius(3));
mask1 = min(mask1,-shapeSphere(g, maskCenter(:,1), maskRadius(1)));

mask = max(mask1,-shapeSphere(g, maskCenter(:,2), maskRadius(2)));

    