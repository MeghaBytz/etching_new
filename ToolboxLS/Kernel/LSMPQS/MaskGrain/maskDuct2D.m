function [mask,g] = maskDuct2D(g)

% function [mask,g] = maskDuct2D(g)
% Create mask (intersection of hyperplanes creates a horizontal duct)
% g - grid structure
    
normal = [ 0.0; -1.0];
point  = [g.min(1);g.min(2) + g.dx(2)];

mask = shapeHyperplane(g,normal,point);

normal = [ 0.0; 1.0];
point  = [g.min(1);g.max(2) - g.dx(2)];

mask = max(mask,shapeHyperplane(g,normal,point));