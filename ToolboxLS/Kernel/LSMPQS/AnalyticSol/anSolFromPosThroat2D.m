function [ data, C ] = anSolFromPosThroat2D(x,g,mask)

% function [ data, C ] = anSolFromPosThroat2D(x,g,mask)
% Analytical solution corresponding to 'x' coordinate for Throat geometry
% whose level set function is 'mask' and grid is 'g' (see maskThroat2D.m 
% for geometry details)
%
% Returns 
% data - level set function that corresponds to the analytical solution
% C - corresponding curvature

%finding effective throat size from data
%central1 = ceil(g.N(1)/2);
%central2 = ceil(g.N(2)/2);
%thr_eff_sz = size( find(mask(central1,:) < 0) );
%thr_eff_sz = thr_eff_sz(1)*thr_eff_sz(2);
%thr_eff_sz = thr_eff_sz*g.dx(2)/2.0;
%R0 = thr_eff_sz;

R = 1.0;
R0 = 0.15; % half of the gap

r1 = 2*R*R*( 1 - sqrt(1 - x*x/(R*R))) - x*x;
r1 = sqrt(r1);

r = R0 + r1;
rc = r*R/(R + R0 - r);
x1 = sqrt(rc*rc - r*r);

if( x < 0 )
    xc = x - x1;
else
    xc = x+x1;
end

C = 1/rc;

Center = [xc; 0.0; 0.0; 0.0 ];

data =  shapeSphere(g, Center, rc);
% this is to get the half-circle
normal = [ 1.0; 0.0];
point = [xc; 0.0];
data = min(data, shapeHyperplane(g,normal,point));

data = max(data,mask);

%plot
figure
lev = [ 0 0 ];
  
[ garbage, hI ] = contour(g.xs{1}, g.xs{2}, data, lev, 'r-');
hold on;
plot([x;x],[0;r],'b-')
plot([xc;x],[0;r],'b-')
[ garbage, hM ] = contourf(g.xs{1}, g.xs{2}, mask, lev, 'k-'); colormap gray;
hs = [ hI; hM];
%set(hs, 'LineWidth',2.0);
plot([x;x],[0;r],'b-')

legend([ hI(1), hM(1)], 'anal.sol.', 'mask');
axis equal
axis(g.axis);

