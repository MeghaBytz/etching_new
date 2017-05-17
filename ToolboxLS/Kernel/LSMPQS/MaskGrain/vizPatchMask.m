
function mask = vizPatchMask(GeomType,GeomType1,PatchX,PatchY,GeomType2,GeomType3)

% function mask = vizPatchMask(GeomType,GeomType1,PatchX,PatchY,GeomType2)
% visualizes specific patchwork of 2 geometries
% GeomType[ - geom. type, see setAvailGeom.m for complete list
% Patchx,patchy       - grid spacing

if( nargin == 4 )
    options = setOptions('GeomType',GeomType,'GeomType1',GeomType1,'PatchX',PatchX,'PatchY',PatchY);
elseif( nargin == 5 )
    options = setOptions('GeomType',GeomType,'GeomType1',GeomType1,'PatchX',PatchX,'PatchY',PatchY,'Geomtype2',GeomType2);
elseif( nargin == 6 )
    options = setOptions('GeomType',GeomType,'GeomType1',GeomType1,'PatchX',PatchX,'PatchY',PatchY,'Geomtype2',GeomType2,'Geomtype3',GeomType3);
else
    error('Incomplete list of input arguments');
end

[AvailGeom, MaskGeom] = setAvailGeom;

[data, data0, g, mask] = initializeDataGrid(options,MaskGeom,0);

if(g.dim == 2)
    figure,
    contourf(g.xs{1}, g.xs{2}, mask, [0 0], 'k-'); colormap gray;
    axis(g.axis); axis image;
    title('Combo');
else error('Plots available in 2D only.');
end