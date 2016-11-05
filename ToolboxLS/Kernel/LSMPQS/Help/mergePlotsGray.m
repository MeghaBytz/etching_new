function mergePlotsGray(d0,d1,d2,d3,g,mask,d1_title,d2_title,d3_title)
% function mergePlotsGray(d0,d1,d2,d3,g,mask,d1_title,d2_title)
% function mergePlotsGray(d0,d1,d2,d3,g,mask,d1_title,d2_title,d3_title)
% d0 data - initial condition
% d1,d2 - mid stages
% d3 final data, if d3_title not provided it will be called 'critical'
% g  - grid
% mask - mask for pore space
% d1_title - title for d1 data
% d2_title - title for d2 data
% Display initial set, mask, mid stages and final set.

%set flip to 1 if you want flip x and y axis on the plot
flip = 1;
if( flip )
    mask = mask';
    d0 = d0'; d1 = d1'; d2 = d2'; d3 = d3';
    tmp = g.xs{1}; g.xs{1} = g.xs{2}; g.xs{2} = tmp;
    g.xs{1} = g.xs{1}';
    g.xs{2} = g.xs{2}';
    tmp1 = g.axis;
    g.axis(1) = g.axis(3);
    g.axis(2) = g.axis(4);
    g.axis(3) = tmp1(1);
    g.axis(4) = tmp1(2);
    tmp2 = g.N(1); g.N(1 ) = g.N(2); g.N(2) = tmp2;
    tmp2 = g.shape(1); g.shape(1) = g.shape(2);  g.shape(2) = tmp2; 
end

  figure;
  lev = [ 0 0 ];
  
  [ garbage, hI ] = contour(g.xs{1}, g.xs{2}, d0, lev, 'k:');
  hold on;
  [ garbage, hM ] = contourf(g.xs{1}, g.xs{2}, mask, lev, 'k-'); colormap gray;
  [ garbage, hD1 ] = contour(g.xs{1}, g.xs{2}, d1, lev, 'Color',[0.6,0.6,0.6]);
  [ garbage, hD2 ] = contour(g.xs{1}, g.xs{2}, d2, lev, 'k:');
  [ garbage, hD3 ] = contour(g.xs{1}, g.xs{2}, d3, lev, 'k-');
  
  
  hs = [ hD1; hD2; hD3];

  set(hs, 'LineWidth',2.0);
   
  if( nargin < 9) d3_title = 'critical'; end
  
  legend([ hI(1), hD1(1), hD2(1), hD3(1)], ...
         {'initial', d1_title, d2_title, d3_title}, 2);
  axis equal
  axis(g.axis);