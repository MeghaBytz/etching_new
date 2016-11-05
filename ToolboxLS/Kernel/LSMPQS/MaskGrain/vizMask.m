
function mask = vizMask(GeomType,segfilename)

% function mask = vizMask(GeomType)
% functionmask -  vizMask(GeomType,segfilename)

% visualizes a specific geometry
% GeomType - geom. type, see setAvailGeom.m for complete list
% segfilename - needs to be provided if GeomType = 'SegFile3DMA'

options = setOptions('GeomType',GeomType);
if( nargin >= 2 )
    %options.doSeal = 1;
    %options.doFlip = 1;
    %options.ICplane = 'y';
    options = setOptions(options,'segFileName',segfilename)
end

[AvailGeom, MaskGeom] = setAvailGeom;

[data, data0, g, mask] = initializeDataGrid(options,MaskGeom,0);

if(g.dim == 2)
    figure,
    contourf(g.xs{1}, g.xs{2}, mask, [0 0], 'k-'); colormap gray;
    axis(g.axis); axis image;
    title(options.GeomType);
elseif( g.dim == 3 )
    % 3D visualization
    x = [1:g.N(1)]; y = [1:g.N(2)]; z = [1:g.N(3)];
    [X Y Z] = meshgrid(x,y,z);
    figure, h = patch(isosurface(X,Y,Z,mask,0));
    set(h,'FaceColor',[0.8 0.8 0.8],'EdgeColor', 'none');
    daspect([1 1 1]); camlight; lighting phong
    view(3)
else
    error('Plots available in 2D/3D only.');
end