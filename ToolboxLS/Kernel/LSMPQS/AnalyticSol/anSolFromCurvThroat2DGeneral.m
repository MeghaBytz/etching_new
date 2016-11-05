function x = anSolFromCurvThroat2DGeneral(C,theta,g,mask,negative,doPlot)

% function [ data, x ] = anSolFromCurvThroat2D(C,theta,g,mask,negative,doPlot)
% Analytical solution corresponding to 'x' coordinate for Throat geometry
% whose level set function is 'mask' and grid is 'g' (see maskThroat2D.m 
% for geometry details)
% theta - contact angle in degrees
%
% Returns 
%  x - position

%finding effective throat size from data
%central1 = ceil(g.N(1)/2);
%central2 = ceil(g.N(2)/2);
%thr_eff_sz = size( find(mask(central1,:) < 0) );
%thr_eff_sz = thr_eff_sz(1)*thr_eff_sz(2);
%thr_eff_sz = thr_eff_sz*g.dx(2)/2.0;
%R0 = thr_eff_sz;


R = 1.0;
R0 = 0.15; % half of the gap

n = numel(C);

if( doPlot )
    figure
    lev = [ 0 0 ];
    title_plot = sprintf('theta %g, C= ',theta);
end
    
for i=1:n
    
 rc = 1/C(i);  % radius of curvature

 g1 = rc*cosd(theta);
 v = rc*sind(theta);
 y = sqrt(R*R - v*v);
 x(i) = (y+g1)*(y+g1) - (R+R0)*(R+R0);


 x(i) = sqrt(x(i));


 if( negative )
    x(i) = -x(i);
 end
 xc = x(i)

 Center = [xc; 0.0; 0.0; 0.0 ]

 data =  shapeSphere(g, Center, rc);
 % this is to get the half-circle
 normal = [ 1.0; 0.0];
 point = [xc; 0.0];
 data = min(data, shapeHyperplane(g,normal,point));

 data = max(data,mask);

 %plot
 if( doPlot )

    [ garbage, hI ] = contour(g.xs{1}, g.xs{2}, data, lev, 'r-');
    hold on;
    
    [ garbage, hM ] = contourf(g.xs{1}, g.xs{2}, mask, lev, 'k-'); colormap gray;
    hold on;
    axis equal
    axis(g.axis);
    
    if (i < n)
        title_plot = sprintf('%s %g, ',title_plot,C(i));
    else
        title_plot = sprintf('%s %g ',title_plot,C(i));
    end    
 end

end
title(title_plot);
hold off
