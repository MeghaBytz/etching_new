function [mask,g] = maskHybrid1(g)

% function [mask,g] = maskHybrid1(g)
% Create mask for a hybrid btw 'Biconic2D' and 'Throat2D'
% g - grid structure
% 
% [-0.5-g.dx(1), 0.5+g.dx(1)] x [-0.3-g.dx(2), 0.3+g.dx(2)]

R = 1.0;
gap = 0.3;

maskCenter = [0; gap/2 + R; 0.0; 0.0 ];
maskRadius = R;
mask = -shapeSphere(g, maskCenter, maskRadius);

normal = [ 0.3; -1.0];
point  = [g.min(1);g.min(2) + g.dx(2)];

mask1 = shapeHyperplane(g,normal,point);

normal = [ -0.3; -1.0];
point  = [g.max(1);g.min(2) + g.dx(2)];

mask1 = min(mask1,shapeHyperplane(g,normal,point));

mask = max(mask,mask1);


    