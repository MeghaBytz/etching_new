function data = createICimb(g,plane,grid_pos)

% data = createICimb(g,plane,grid_pos)
% Create initial conditions for imbibition
% Nonwetting fluid is filling most of the space.
% g - grid
% plane - 'x' (plane x=0), 'y' (plane y=0) or 'z' (plane z = 0)
% where - integer identifying at which grid point shall we put the plane
%         anywhere from 1 to ( g.N(i)-1 ) (i=1,2,3 depending on the plane)

if( grid_pos < 1) grid_pos = 1; end
    
if( g.dim == 2 )
    % in 2D hyperplane x=0 at (g.min+g.dx,0,0).
    if( plane == 'x')
        normal = [ -1.0; 0.0];
        point  = [g.min(1) + grid_pos*g.dx(1);g.min(2)];
    else
        normal = [0.0;-1.0];
        point  = [g.min(1);g.min(2) + g.dx(2)];
    end     
else
    if( plane == 'z')
        normal = [0.0; 0.0; -1.0];
        point  = [0.0; 0.0;g.min(3) + grid_pos*g.dx(3)];
    elseif( plane == 'y')
        normal = [0.0; -1.0;0.0];
        point  = [0.0;g.min(2) + grid_pos*g.dx(2);0.0];
    else
        normal = [-1.0; 0.0; 0.0];
        point  = [g.min(1) + grid_pos*g.dx(1);0.0;0.0];
    end
end
data = shapeHyperplane(g,normal,point); 