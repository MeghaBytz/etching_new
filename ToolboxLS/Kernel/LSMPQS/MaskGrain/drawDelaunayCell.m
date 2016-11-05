function drawDelaunayCell(dxIn) 
%---------------------------------------------------------------------------
% Create mask (void space btw 4 spheres and 3 hyperplanes)
%---------------------------------------------------------------------------
r = 1;  % sphere radius
n = 4;  % number of masking spheres

% 4 spheres of radius r have centers at columns of maskCenter
eps = 0.04;
a = 2*r;
v = a*sqrt(3)/2.0;
k = a *sqrt(2.0/3.0);
maskCenter = [0 0 0; (a+eps) 0 0; a/2 (v+eps) 0; a/2 v/3 (k+eps)];
maskCenter = maskCenter';


g.dim = 3;

if( nargin  == 0)
    dxIn = 0.02;
end
nx = ceil(1/dxIn); % need to avoid rounding off errors
g.dx = 1/nx; %input value

g.min = 0 - 2*g.dx;
g.max = 2;
g.bdry = @addGhostExtrapolate;
g = processGrid(g)
 
% let mask be negative inside the 4 spheres
mask = shapeSphere(g, maskCenter(:,1), r);
for i=2:n
  mask = min(mask,shapeSphere(g, maskCenter(:,i),r));  
end



A = maskCenter(:,1)
B = maskCenter(:,2)
C = maskCenter(:,3)
D = maskCenter(:,4)

normal = cross(B-A,D-A);
planeABD = shapeHyperplane(g,normal,A);

normal = cross(D-A,C-A);
planeADC = shapeHyperplane(g,normal,A);

normal = cross(C-A,B-A);
planeACB = shapeHyperplane(g,normal,A);

normal = cross(C-B,D-B);
planeBCD = shapeHyperplane(g,normal,B);

tetra = max(planeABD,planeADC);  
tetra = max(tetra,planeACB);  
tetra = max(tetra,planeBCD);  
 
%intersect spheres and tetrahedron
mask = max(mask,tetra);


% 3D visualization
x = [1:g.N(1)]; y = [1:g.N(2)]; z = [1:g.N(3)];
[X Y Z] = meshgrid(x,y,z);
figure, h = patch(isosurface(X,Y,Z,mask,0));
set(h,'FaceColor',[0.9 0.9 0.9],'EdgeColor', 'none');
daspect([1 1 1]); camlight; axis off;
view(3)
