function testGradient(data,grid,mask)

% function testGradient(data,grid,mask)
% data - level set function (in 2d)
% grid - grid structure
% mask - mask level set function
% Plots pointwise magnitude of gradient for data and mask 

 phi_x = centeredFirstSecond(grid, data, 1);
 phi_y = centeredFirstSecond(grid, data, 2);
 mask_x = centeredFirstSecond(grid, mask, 1);
 mask_y = centeredFirstSecond(grid, mask, 2);

 prod1 = phi_x .* phi_x + phi_y .* phi_y;
 prod1 = sqrt(prod1);
 figure, plot3(grid.xs{1}, grid.xs{2},prod1); title('Gradient magnitude for data');
 
 prod2 = mask_x .* mask_x + mask_y .* mask_y;
 prod2 = sqrt(prod2);
 figure, plot3(grid.xs{1}, grid.xs{2},prod2); title('Gradient magnitude for mask');
 
 prod3 = phi_x .* mask_x + phi_y .* mask_y;
 angle_term = prod3 ./ (prod1 .* prod2);
 
 problem_pts = find( prod1 == 0);
 angle_term(problem_pts) = 0;
 problem_pts = find( prod2 == 0);
 angle_term(problem_pts) = 0;
 
 figure, plot3(grid.xs{1}, grid.xs{2},acos(angle_term)); title('Angle between gradients');