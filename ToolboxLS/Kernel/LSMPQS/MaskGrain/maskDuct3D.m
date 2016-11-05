function [mask,g] = MaskDuct3D(g)
%--------------------------------------------------------------------------
% Create mask (intersection of hyperplanes creates a horizontal duct)
%--------------------------------------------------------------------------

    % y-direction duct sides
    normal = [ 0.0; -1.0; 0.0];
    point  = [g.min(1);g.min(2) + g.dx(2); g.min(3)];
   
    mask = shapeHyperplane(g,normal,point);

    normal = [ 0.0; 1.0;0.0];
    point  = [g.min(1);g.max(2) - g.dx(2); g.min(3)];

    mask = max(mask,shapeHyperplane(g,normal,point));
    
    % z-direction duct sides
    normal = [ 0.0; 0.0; -1.0];
    point  = [g.min(1);g.min(2); g.min(3) + g.dx(3)];
   
    mask = max(mask,shapeHyperplane(g,normal,point));

    
    normal = [ 0.0; 0.0; 1.0];
    point  = [g.min(1);g.min(2); g.max(3) - g.dx(3)];

    mask = max(mask,shapeHyperplane(g,normal,point));