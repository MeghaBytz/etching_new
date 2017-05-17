function  [mask,g] = maskSquarePore(g)

% function  [mask,g] = maskSquarePore(g)
% Create square nodal pore typical for network algorithms.

width = 0.35 + 2*g.dx(1); %leaves gap og 0.3 in the middle
normal = [ 0.0; -1.0];
point  = [g.min(1);g.min(2) + width];

mask = shapeHyperplane(g,normal,point);

normal = [ 0.0; 1.0];
point  = [g.min(1);g.max(2) - width];

mask = max(mask,shapeHyperplane(g,normal,point));


normal = [-1.0; 0];
point  = [g.min(1) + width;g.min(2)];

mask1 = shapeHyperplane(g,normal,point);

normal = [ 1.0; 0.0];
point  = [g.max(1) - width; g.min(2)];

mask1 = max(mask1,shapeHyperplane(g,normal,point));

mask = min(mask,mask1);
   
