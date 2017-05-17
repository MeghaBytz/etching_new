function [angle median_angle mean_angle] = explore_contact_angle_vert_main(data1,data2,nx,ny,nz,isolevel,do_vis,color)

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
v1_sz = size(V1); v1_sz = v1_sz(1)
f1_sz = size(F1); f1_sz = f1_sz(1);

% fluid 2 boundary
[ F2 V2 ] = isosurface(g.xs{1},g.xs{2},g.xs{3},data2,isolevel);
N2 = isonormals(data2,V2);
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

%find vertices on the fluid1-fluid2 interface boundary - that's where we
%want to evaluate contact angle
V = find_bdry_vertices(V1,F_FLD_INTFC);
% vertices interior to fluid1-fluid2 interface
Vint = setdiff(V_FLD_INTFC,V,'rows');

mask = min(data1,data2);  %pore space level set function
mask = -mask;             %grain space level set function
[ FM VM ] = isosurface(g.xs{1},g.xs{2},g.xs{3},mask,isolevel);

%puzzle - if I use N1 = isonormals(g.xs{1},g.xs{2},g.xs{3},data1,V); below
%I get NaNs in my normal values?!
N1 = isonormals(data1,V);
N2 = isonormals(data2,V);
NM = isonormals(mask,V);

fprintf('Angle btw the first fluid normals and grain normals at contact points');
angle = calculate_angle_vec(N1,NM,do_vis);
median_angle = median(angle);
mean_angle = mean(angle);
fprintf('Angle btw the second fluid normals and grain normals at contact points');
calculate_angle_vec(N2,NM,do_vis);

fprintf('Angle btw the first fluid normals and second fluid normals at contact pts- expected close to 180');
calculate_angle_vec(N1,N2,do_vis);

N1int = isonormals(data1,Vint);
N2int = isonormals(data2,Vint);
fprintf('Angle btw the first fluid normals and second fluid normals at interior fluid intfc pts');
calculate_angle_vec(N1int,N2int,do_vis);



