function [mask,g] = rhomboidSpheres3D()
r = 1;
phi = 90;
dx = 0.06;      %seems dx = 0.1 doesn't give a good shape, especially for lower angles

cent_x = [0  2*r 2*r*cosd(phi) 2*r+2*r*cosd(phi)];
cent_y = [0 0 2*r*sind(phi) 2*r*sind(phi)];
cent_z = [0 0 0 0];

x_min = min(cent_x)-r;
x_max = max(cent_x)+r;
y_min = floor(min(cent_y)-r);
y_max = ceil(max(cent_y)+r);
z_min = min(cent_z)-r;
z_max = max(cent_z)+r;

[X Y Z] = meshgrid(x_min:dx:x_max,y_min:dx:y_max, z_min:dx:z_max);

d1 = ((X-cent_x(1)).^2+(Y-cent_y(1)).^2+(Z-cent_z(1)).^2).^0.5 - r;
d2 = ((X-cent_x(2)).^2+(Y-cent_y(2)).^2+(Z-cent_z(2)).^2).^0.5 - r;
d3 = ((X-cent_x(3)).^2+(Y-cent_y(3)).^2+(Z-cent_z(3)).^2).^0.5 - r;
d4 = ((X-cent_x(4)).^2+(Y-cent_y(4)).^2+(Z-cent_z(4)).^2).^0.5 - r;

mask = -1*min(min(d1,d2),min(d3,d4));
mask((Y>X*tand(phi))&(X>cent_x(1))&(X<cent_x(3)))=1;
mask((Y<X*tand(phi))&(X>cent_x(2))&(X<cent_x(4)))=1;
mask(X>=max(cent_x)) = 1;
mask(Y>=max(cent_y)) = 1;
mask(X<=min(cent_x)) = 1;
mask(Y<=min(cent_y)) = 1;
%mask(Z>=z_max-r) = -1;
%mask(Z<=z_min+r) = -1;

g.dim = 3;
g.min = [y_min; x_min; z_min];
g.max = [y_max; x_max; z_max];
g.dx = dx;
g.bdry = @addGhostExtrapolate;
g = processGrid(g);

%visualizeLevelSet(g,mask,'surface',0);
%axis equal

end