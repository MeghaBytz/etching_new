function [angle median_angle mean_angle] = explore_contact_angle_main(data1,data2,nx,ny,nz,isolevel,do_vis,color)

% [angle median_angle mean_angle] = explore_contact_angle_main(data1,data2,nx,ny,nz,isolevel,do_vis,color)
% data1,2 - [nx ny nz] arrays containing fluid1,2 data
%              the first file is treated as the main fluid of interest
% do_vis - set to 1 if you want to visualaize result
% isolevel - indicates level set of the interface ni data1 and data2
% angle  = contact angle calculated at vertices of contact btw
%          fluid1-fluid2 isosurface and fluid1-grain isosurface
% median_angle,mean_angle - median and mean of the array angle


% create grid, needed for curvature computation (with LS Toolbox function)
g.dim = 3;
g.min = 1;
g.dx = 1; % irrelevant for now
g.max = g.dx*[nx; ny; nz];
g.bdry = @addGhostExtrapolate;
g = processGrid(g)

% fluid 1 boundary
[ F1 V1 ] = isosurface(g.xs{1},g.xs{2},g.xs{3},data1,isolevel);
%N1 = isonormals(g.xs{1},g.xs{2},g.xs{3},data1,V1)
v1_sz = size(V1); v1_sz = v1_sz(1)
f1_sz = size(F1); f1_sz = f1_sz(1);

% fluid 2 boundary
[ F2 V2 ] = isosurface(g.xs{1},g.xs{2},g.xs{3},data2,isolevel);
v2_sz = size(V2); v2_sz = v2_sz(1)
f2_sz = size(F2); f2_sz = f2_sz(1);

%find vertices and triangles on fluid1-fluid2 interface
[ V_FLD_INTFC i1 i2] = intersect(V1,V2,'rows');

%find triangles that have all 3 vertices in V_FLD_INTFC
v1_fld = zeros(v1_sz,1); v1_fld(i1) = 1;
f_fld = find((v1_fld(F1(:,1)) & v1_fld(F1(:,2)) & v1_fld(F1(:,3))));
F_FLD_INTFC = F1(f_fld,:);

%isolate all triangles that are not in F_FLD_INTFC - for visualization purpose
f_grain = find( (v1_fld(F1(:,1))==0) | (v1_fld(F1(:,2))==0) | (v1_fld(F1(:,3))==0) );
F_GRAIN = F1(f_grain,:);

% visualize
if( do_vis )
    figure, h_fld = patch('Faces',F_FLD_INTFC,'Vertices',V1);
    set(h_fld,'FaceColor',color,'EdgeColor', 'none');
    daspect([1 1 1]); camlight; lighting phong;
    view(3), hold on
    
    fprintf('\nFluid1-grain boundary is shown in gray, fluid1-fluid2 surface is red or green.');

    h_grain = patch('Faces',F_GRAIN,'Vertices',V1);
    size(F_GRAIN)
    set(h_grain, 'FaceColor', [0.8 0.8 0.8] ,'EdgeColor', 'none'); % gray color
    daspect([1 1 1]); camlight; lighting phong;
    view(3), hold off
end

[ADJ_FLD ADJ_GRAIN ] = find_adjacent_faces(F_FLD_INTFC,F_GRAIN);
angle = calculate_face_angles(V1,F_FLD_INTFC,ADJ_FLD,F_GRAIN,ADJ_GRAIN);

if( do_vis ) figure, hist(angle,100);  end
num_of_bdry_points = size(angle);
median_angle  = median(angle)
mean_angle = mean(angle)

% visualize adjacent faces
if( do_vis )    
    figure, h_fld = patch('Faces',F_FLD_INTFC(ADJ_FLD,:),'Vertices',V1);
    set(h_fld,'FaceColor',color,'EdgeColor', 'none');
    daspect([1 1 1]); camlight; lighting phong;
    view(3), hold on
    
    
    h_grain = patch('Faces',F_GRAIN(ADJ_GRAIN,:),'Vertices',V1);
    size(F_GRAIN)
    set(h_grain, 'FaceColor', [0.8 0.8 0.8] ,'EdgeColor', 'none'); % gray color
    daspect([1 1 1]); camlight; lighting phong;
    view(3), hold off
end


