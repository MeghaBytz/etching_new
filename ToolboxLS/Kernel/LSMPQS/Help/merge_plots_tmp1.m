function mergePlots(d0,d1,g,mask)
% d0 data - initial condition
% d1,d2,d3 data - three different stages
% g  - grid
% mask - mask for pore space
% Display initial set, mask, minimum over time, and final set.
  figure;
  lev = [ 0 0 ];
  
  [ garbage, hI ] = contour(g.xs{1}, g.xs{2}, d0, lev, 'r-');
  hold on;
  [ garbage, hM ] = contourf(g.xs{1}, g.xs{2}, mask, lev, 'k-'); colormap gray;
  [ garbage, hD1 ] = contour(g.xs{1}, g.xs{2}, d1, lev, 'k:');
  %[ garbage, hD2 ] = contour(g.xs{1}, g.xs{2}, d2, lev, 'm-');
  %[ garbage, hD3 ] = contour(g.xs{1}, g.xs{2}, d3, lev, 'k:');
  
  
  hs = [ hI; hD1; hM];

  set(hs, 'LineWidth',2.0);

  legend([ hI(1), hD1(1)], 'data-step4', 'data-a0+da');
  axis equal
  axis(g.axis);