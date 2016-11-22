function [data, mask, mask_vel, init_level] = initializeGeom(g,options)

%---------------------------------------------------------------------------
% Create initial conditions (a rectangle).
widths = (g.max-g.min)+2*g.dx;
center = g.min + 0.5*(g.max - g.min);
center(end) = center(end) - 0.2;
data = shapeRectangleByCenter(g,center,widths);
init_level = center(end)+widths(end)/2;

%---------------------------------------------------------------------------
% Create mask 

if(options.doMask)
    half_width =0.1;
    thickness = 6*g.dx;
    mask_widths = [(g.max(1:end-1) - g.min(1:end-1) + 2*g.dx(1)); thickness];
    mask_centers = [(g.min(1:end-1) + 0.5*(g.max(1:end-1) - g.min(1:end-1))); init_level];

    mask = shapeRectangleByCenter(g,mask_centers,mask_widths);
    mask_vel = shapeRectangleByCenter(g,mask_centers,[mask_widths(1);mask_widths(2) + 2*g.dx]);

    hole_centers = [(g.min(1:end-1) + 0.5*(g.max(1:end-1) - g.min(1:end-1))); init_level];
    hole_widths = [half_width*2*(1:g.dim-1); thickness+2*g.dx];
    hole = shapeRectangleByCenter(g,hole_centers,hole_widths);
    hole_vel = shapeRectangleByCenter(g,hole_centers,[hole_widths(1); hole_widths(2) + 2*g.dx]);

    mask = shapeDifference(mask,hole);
    mask_vel = shapeDifference(mask_vel, hole_vel);

    % we need the complement of the masked region.
    mask = -mask;
    mask_vel = -mask_vel;

    % Need to ensure that the initial conditions satisfy the mask.
    data_new = max(data, mask);

    data = data_new;
else
    mask = []; mask_vel = [];
end
