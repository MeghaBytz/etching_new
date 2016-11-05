function [mask,g] = maskBiconic2D(g)

% function [mask,g] = maskBiconic2D(g)
% Create mask for a 2D biconical section
% g - grid structure
% This biconical section is specific for geometry
% [-0.5-g.dx(1), 0.5+g.dx(1)] x [-0.3-g.dx(2), 0.3+g.dx(2)]

normal = [ 0.3; -1.0];
point  = [g.min(1);g.min(2) + g.dx(2)];

mask = shapeHyperplane(g,normal,point);

normal = [ -0.3; -1.0];
point  = [g.max(1);g.min(2) + g.dx(2)];

mask = min(mask,shapeHyperplane(g,normal,point));

normal = [ 0.3; 1.0];
point  = [g.min(1);g.max(2) - g.dx(2)];

mask1 = shapeHyperplane(g,normal,point);

normal = [ -0.3; 1.0];
point  = [g.max(1);g.max(2) - g.dx(2)];

mask1 = min(mask1,shapeHyperplane(g,normal,point));

mask = max(mask,mask1);
    
    