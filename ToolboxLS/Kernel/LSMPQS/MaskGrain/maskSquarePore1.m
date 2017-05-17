function  [mask,g] = maskSquarePore1(g)

% function  [mask,g] = maskSquarePore(g)
% Create square nodal pore typical for network algorithms.

width1 = 0.35 + 2*g.dx(1); %leaves gap og 0.3 in the middle
width = 0.3 + 2*g.dx(1); %leaves gap og 0.4 in the middle

normal = [ 0.0; -1.0];
point  = [g.min(1);g.min(2) + width];

mask = shapeHyperplane(g,normal,point);

normal = [ 0.0; 1.0];
point  = [g.min(1);g.max(2) - width];

mask = max(mask,shapeHyperplane(g,normal,point));


normal = [-1.0; 0];
point  = [g.min(1) + width1;g.min(2)];

mask1 = shapeHyperplane(g,normal,point);

normal = [ 1.0; 0.0];
point  = [g.max(1) - width1; g.min(2)];

mask1 = max(mask1,shapeHyperplane(g,normal,point));

mask = min(mask,mask1);
   
