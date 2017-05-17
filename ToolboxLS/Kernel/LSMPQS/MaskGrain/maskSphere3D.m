function [mask,g] = maskSphere3D(g)    
%---------------------------------------------------------------------------
% Create mask (void space btw 4 spheres and 3 hyperplanes)
%---------------------------------------------------------------------------
r = 1;  % sphere radius
n = 4;  % number of masking spheres

% 4 spheres of radius r have centers at columns of maskCenter
a = 2*r;
v = a*sqrt(3)/2.0;
k = a *sqrt(2.0/3.0);
maskCenter = [0 0 0; a 0 0; a/2 v 0; a/2 v/3 k];
maskCenter = maskCenter';

% The moving set can be anywhere outside the masked region.
%  '-' sign is to get the complement of the masked region.
mask = - shapeSphere(g, maskCenter(:,1), r);


for i=2:n
  mask = max(mask,-shapeSphere(g, maskCenter(:,i),r));  
end
figure, visualizeLevelSet(g, mask, 'surface', 0);

point0 = maskCenter(:,4);
for i=1:(n-1)
   point1 = maskCenter(:,i);
   if ( i == 3) point2 = maskCenter(:,1);
   else point2 = maskCenter(:,i+1);
   end
   
   A = point1 - point0;
   B = point2 - point0;
   normal = cross(A,B);
  normal = normal/norm(normal);
   
  mask = max(mask,shapeHyperplane(g,normal,point0));
end
