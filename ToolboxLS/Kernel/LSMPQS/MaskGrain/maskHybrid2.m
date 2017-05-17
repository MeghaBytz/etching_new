function [mask,g] = maskHybrid2(g)

% function [mask,g] = maskHybrid2(g)
% Create mask for a hybrid - piecewise biconical and duct.
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
normal = [ 0; -1.0];
point  = [-0.3;-0.15];

mask3 = shapeHyperplane(g,normal,point);

normal = [ 0; 1.0];
point  = [-0.4;0.15];

mask3 = max(mask3,shapeHyperplane(g,normal,point));

normal = [ -1.0; 0];
point  = [-0.1-g.dx(1)/2.0;0];

mask3 = max(mask3,shapeHyperplane(g,normal,point));

%final
mask = min(mask1,mask2);
mask = min(mask,mask3);


    