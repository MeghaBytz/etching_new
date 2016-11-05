function data = rectangleIC(g)

% data = rectangleIC(g)
% g - grid

diff = g.max - g.min;
p = g.min + diff/3.0;
q = g.min + 2.0*diff/3.0;


if( g.dim == 2 )
    % in 2D hyperplane x=0 at (g.min+g.dx,0,0).
    normal = [ -1.0; 0.0];
    point1  = [p(1);p(2)];
    data = shapeHyperplane(g,normal,point1); 
    
    normal = - normal;
    point2  = [q(1);q(2)];
    data = max(data,shapeHyperplane(g,normal,point2)); 
    
    normal = [0;-1.0];
    %data = max(data,shapeHyperplane(g,normal,point1)); 
    
    normal = -normal;
    %data = max(data,shapeHyperplane(g,normal,point2)); 
else
    disp('Fix rectangleIC.m for 3D');
end