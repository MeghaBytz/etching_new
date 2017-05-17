function plotGradient2D(data,grid,mask)

% function plotGradient2D(data,grid,mask)
% plots (interpolated) gradient values for level set function data at the
% zero level set

 phi_x = centeredFirstSecond(grid, data, 1);
 phi_y = centeredFirstSecond(grid, data, 2);
 
 mask_x = centeredFirstSecond(grid, mask, 1);
 mask_y = centeredFirstSecond(grid, mask, 2);

 %figure, quiver(grid.xs{1}, grid.xs{2},phi_x,phi_y);
 
 C = contourc(data,[0 0]); % find a 0-level set and output coords to matrix C
 C = getContourPoints(C); C = C'; %isolate coords of points on the level set
 x = C(:,1);
 y = C(:,2);
 phi_x_interp = interp2(phi_x,x,y);
 phi_y_interp = interp2(phi_y,x,y);
 
 mask_x_interp = interp2(mask_x,x,y);
 mask_y_interp = interp2(mask_y,x,y);
 
 translate = ones(size(x));
 x = grid.min(2)*translate + (x - translate )*grid.dx(2);
 y = grid.min(1)*translate + (y - translate )*grid.dx(1);
 
 figure
 contourf(grid.xs{1}, grid.xs{2}, mask, [0 0] , 'k-');colormap gray; hold on;
 quiver(y,x,phi_x_interp,phi_y_interp); hold on;
 quiver(y,x,mask_x_interp,mask_y_interp,'m');
 contour(grid.xs{1}, grid.xs{2}, data, [0 0] , 'r-');
 axis(grid.axis); axis image;
 hold off
