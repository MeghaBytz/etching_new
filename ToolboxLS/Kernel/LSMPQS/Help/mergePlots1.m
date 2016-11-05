function mergePlots1(d0,d1,d2,d3,d4,d5,g,mask,d1_title,d2_title,d3_title,d4_title,d5_title)
% function mergePlots1(d0,d1,d2,d3,d4,d5,g,mask,d1_title,d2_title,d3_title)
% function mergePlots1(d0,d1,d2,d3,d4,d5,g,mask,d1_title,d2_title,d3_title,d4_title)
% function mergePlots1(d0,d1,d2,d3,d4,d5,g,mask,d1_title,d2_title,d3_title,d4_title,d5_title)
% d0 data - initial condition
% d1,d2,d3-  data - three different stages
% d4 - critical stage; if d4_title not provided it will be called
% 'critical'
% d5 - final stage; if d5_title not provided it will be called
% 'final'
% g  - grid
% mask - mask for pore space
% Display initial set, mask, minimum over time, and final set.
  figure;
  lev = [ 0 0 ];
  
  [ garbage, hI ] = contour(g.xs{1}, g.xs{2}, d0, lev, 'b-');
  hold on;
  [ garbage, hM ] = contourf(g.xs{1}, g.xs{2}, mask, lev, 'k-'); colormap gray;
  [ garbage, hD1 ] = contour(g.xs{1}, g.xs{2}, d1, lev, 'r-');
  [ garbage, hD2 ] = contour(g.xs{1}, g.xs{2}, d2, lev, 'g-');
  [ garbage, hD3 ] = contour(g.xs{1}, g.xs{2}, d3, lev, 'c-');
  [ garbage, hD4 ] = contour(g.xs{1}, g.xs{2}, d4, lev, 'm-');
  [ garbage, hD5 ] = contour(g.xs{1}, g.xs{2}, d5, lev, 'Color',[1.0 0.5 0]);
 
  
  hs = [ hI; hD1; hD2; hD3; hD4; hD5; hM];

  set(hs, 'LineWidth',2.0);

  if( nargin < 12) d4_title = 'critical'; end
  if( nargin < 13) d5_title = 'final'; end
  
  legend([ hI(1), hD1(1), hD2(1), hD3(1),hD4(1),hD5(1)], ...
         {'initial', d1_title, d2_title, d3_title,d4_title,d5_title}, 2);
  axis equal
  axis(g.axis);