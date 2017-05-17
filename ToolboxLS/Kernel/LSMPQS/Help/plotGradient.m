function plotGradient(data,grid)

% function plotGradient(data,grid)
% plot pointwise level set function 'data' gradient vectors

 phi_x = centeredFirstSecond(grid, data, 1);
 phi_y = centeredFirstSecond(grid, data, 2);
 
 mask_x = centeredFirstSecond(grid, mask, 1);
 mask_y = centeredFirstSecond(grid, mask, 2);

 figure, quiver(grid.xs{1}, grid.xs{2},phi_x,phi_y);
 