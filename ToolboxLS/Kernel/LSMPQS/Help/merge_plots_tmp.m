function mergePlots(d0,d1,d2,d3,d4,g,mask)
% d0 data - initial condition
% d1,d2,d3 data - three different stages
% g  - grid
% mask - mask for pore space
% Display initial set, mask, minimum over time, and final set.
  figure;
  lev = [ 0 0 ];
  
  % garbage, hI ] = contour(g.xs{1}, g.xs{2}, dinit, lev, 'b-');
  [ garbage, hD0 ] = contour(g.xs{1}, g.xs{2}, d0, lev, 'b');
  hold on;
  [ garbage, hM ] = contourf(g.xs{1}, g.xs{2}, mask, lev, 'k-'); colormap gray;
  [ garbage, hD1 ] = contour(g.xs{1}, g.xs{2}, d1, lev, 'r-');
  [ garbage, hD2 ] = contour(g.xs{1}, g.xs{2}, d2, lev, 'm-');
  [ garbage, hD3 ] = contour(g.xs{1}, g.xs{2}, d3, lev, 'g-');
  [ garbage, hD4 ] = contour(g.xs{1}, g.xs{2}, d4, lev, 'y-');
  
  
  hs = [ hD0;hM;hD1; hD2; hD3];

  set(hs, 'LineWidth',2.0);

  legend([ hD0(1), hD1(1),hD2(1),hD3(1),hD4(1)], 'init','step1', 'step63','step64','step65');
  axis equal
  axis(g.axis);