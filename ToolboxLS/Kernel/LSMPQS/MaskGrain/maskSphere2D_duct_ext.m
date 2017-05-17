function  [mask,g] = maskSphere2D_duct_ext(g)

% function  [mask,g] = maskSphere2D_new(g)
% Create mask that corresponds to void space btw 3 equal spheres.
% Spheres of radius 1.0 aee created (so that calculations are
% automatically "normalized" by sphere radius)

%also creating duct extentions at throats.

N = 3;  % number of masking spheres

maskCenter(:,1) = [0; 1.7];
maskRadius(1) = 1.0;

maskCenter(:,2) = [0; -0.6];
maskRadius(2) = 1.0;

maskCenter(:,3) = [1.75; 0.5];
maskRadius(3) = 1.0;

% The moving set can be anywhere outside the masked region.
%  '-' sign is to get the complement of the masked region.
mask = - shapeSphere(g, maskCenter(:,1), maskRadius(1));
for i=2:N
   mask = max(mask,-shapeSphere(g, maskCenter(:,i), maskRadius(i)));
end

% add duct extentions to throats
next(1) = 2; next(2) = 3; next(3) = 1;

for(i=1:3)
    j = next(i);
    v(:,i) = maskCenter(:,j) - maskCenter(:,i); 
    normv(i) = norm( v(:,i) );
    n(:,i) = [v(2,i); -v(1,i) ];
    Pi = maskCenter(:,i) + (maskRadius(i)/normv(i))*v(:,i);
    Pj = maskCenter(:,j) - (maskRadius(j)/normv(i))*v(:,i);

    mask1 = min( shapeHyperplane(g,n(:,i),Pi), shapeHyperplane(g,-v(:,i),Pi));
    mask2 = min( shapeHyperplane(g,n(:,i),Pj), shapeHyperplane(g, v(:,i),Pj));
    mask = max(mask,mask1);
    mask = max(mask,mask2);
end
%mask = shapeHyperplane(g,-v(:,1),P1);