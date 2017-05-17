function [ data, x ] = anSolFromCurvBiconic2D(C,g,mask,negative,doPlot)

% function [ data, x ] = anSolFromCurvBiconic2D(C,g,mask)
% Analytical solution corresponding to 'x' coordinate for Biconical geometry
% whose level set function is 'mask' and grid is 'g' (see maskBiconic2D.m 
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

R0 = 0.15;
R1 = 0.3; 
R2 = R1;
L1 = 0.5+g.dx(1);
L2 = L1;

tan_theta1 = -(R0-R1)/L1;
cos_theta1 = 1/sqrt(1+tan_theta1*tan_theta1);
sin_theta1 = tan_theta1 * cos_theta1;

tan_theta2 = -(R0-R2)/L2;
cos_theta2 = 1/sqrt(1+tan_theta2*tan_theta2);
sin_theta2 = tan_theta2 * cos_theta2;

rc = 1/C;

if( negative )
   r = rc*cos_theta1;
   x = (R0-r)/tan_theta1;
   xc = x - rc*sin_theta1;
else
    r = rc*cos_theta2;
    x = (r-R0)/tan_theta2;
    xc = x + rc*sin_theta1;
end

Center = [xc; 0.0; 0.0; 0.0 ];

data =  shapeSphere(g, Center, rc);
% this is to get the half-circle
normal = [ 1.0; 0.0];
point = [xc; 0.0];
data = min(data, shapeHyperplane(g,normal,point));

data = max(data,mask);

if( doPlot )
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
end
