function [mask,g] = maskHybrid3(g)

% function [mask,g] = maskHybrid3(g)
% Create mask for a hybrid - piecewise biconical and throat.
% g - grid structure
% 
% [-0.5-g.dx(1), 0.3+g.dx(1)] x [-0.4-g.dx(2), 0.4+g.dx(2)]


% part 1
normal = [ 0.3; -0.3];
point =  [-0.5;-0.4];
mask1 = shapeHyperplane(g,normal,point);

normal = [ 0.3; 0.3];
point  = [-0.5;0.4];

mask1 = max(mask1,shapeHyperplane(g,normal,point));

normal = [ 1.0; 0];
point  = [-0.3+g.dx(1)/2.0;0];

mask1 = max(mask1,shapeHyperplane(g,normal,point));

%part 2
normal = [ 0.3; -1.2];
point  = [-0.3;-0.2];

mask2 = shapeHyperplane(g,normal,point);

normal = [ 0.3; 1.2];
point  = [-0.3; 0.2];

mask2 = max(mask2,shapeHyperplane(g,normal,point));

normal = [ -1.0; 0];
point  = [-0.3;0];


mask2 = max(mask2,shapeHyperplane(g,normal,point));

normal = [ 1.0; 0];
point  = [-0.1;0];

mask2 = max(mask2,shapeHyperplane(g,normal,point));

%part 3

R = 1.0;
gap = 0.3;

maskCenter(:,1) = [-0.1; -gap/2 - R; 0.0; 0.0 ];
maskRadius(1) = R;

maskCenter(:,2) = [-0.1; gap/2 + R; 0.0; 0.0 ];
maskRadius(2) = R;

mask3 = - shapeSphere(g, maskCenter(:,1), maskRadius(1));
mask3 = max(mask3,-shapeSphere(g, maskCenter(:,2), maskRadius(2)));

normal = [ -1.0; 0];
point  = [-0.1-g.dx(1)/2.0;0];

mask3 = max(mask3,shapeHyperplane(g,normal,point));

%final
mask = min(mask1,mask2);
mask = min(mask,mask3);


    