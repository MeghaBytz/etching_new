
function vizMaskIC(GeomType,ICplane,segfilename)

% function vizMaskIC(GeomType,ICplane)
% function vizMaskIC(GeomType,ICplane,segfilename)
% visualizes a specific geometry
% GeomType - geom. type, see setAvailGeom.m for complete list
% ICplane  - initial condition to be visualized
% segfilename - name of the segmented file for geomtype 'Segfile3dma'

options = setOptions('GeomType',GeomType,'Icplane',ICplane);
if( nargin >= 3 )
    options = setOptions(options,'segFileName',segfilename,'doSeal',1);
    disp('doSeal options set to 1');
end
[AvailGeom, MaskGeom] = setAvailGeom;

[data, data0, g, mask] = initializeDataGrid(options,MaskGeom,0);

if(g.dim == 2)
    figure,
    contourf(g.xs{1}, g.xs{2}, mask, [0 0], 'k-'); colormap gray;
    hold on;
    contour(g.xs{1}, g.xs{2}, data0, [0 0], 'b-');
    axis(g.axis); axis image;
    title(options.GeomType);
else error('Plots available in 2D only.');
end