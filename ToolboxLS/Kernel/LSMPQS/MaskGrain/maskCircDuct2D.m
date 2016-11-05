function  [mask,g] = maskCircDuct2D(g)

% function  [mask,g] = maskCircDuct2D(g)

R = 1.0;
gap = 0.3;

maskCenter(:,1) = [0; -gap/2 + R; 0.0; 0.0 ];
maskRadius(1) = R;

maskCenter(:,2) = [0; gap/2 + R; 0.0; 0.0 ];
maskRadius(2) = R;


mask = shapeSphere(g, maskCenter(:,1), maskRadius(1));
mask = max(mask,-shapeSphere(g, maskCenter(:,2), maskRadius(2)));
    