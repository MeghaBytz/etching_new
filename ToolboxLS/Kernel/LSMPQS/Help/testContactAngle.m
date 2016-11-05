function testContactAngle

% function testContactAngle

run('~/ToolboxLS-1.0/Examples/addPathToKernel.m');

do_vis = 1;

%---------Set the geometry---------------------
dxIn = 0.08;

g.dim = 3;
g.min = 0 - 2*dxIn;
nx = ceil(1/dxIn); % need to avoid rounding off errors
g.dx = 1/nx; %input value
g.max = +1 + 2*dxIn; % extra layer is added to each volume side
g.bdry = @addGhostExtrapolate;
g = processGrid(g)

%[mask,g] = MaskDuct3D(g);
Center = [g.min(1)+g.dx(1); 0.5; 0.5; 0.0 ];
Radius = 0.5;

mask = shapeCylinder(g,[1.0;0.0; 0.0],[0.5; 0.5; 0.5; 0.0 ],Radius);

% fluids will have spherical interface
% fluid1 is outside the sphere - wetting fluid
data1 = - shapeSphere(g, Center, Radius);
data1 = max(mask,data1);

% fluid2 is inside the spheer
data2 = shapeSphere(g, Center, Radius);
data2 = max(mask,data2);

%-------Compute interphase surfaces-----------------
isolevel = 0;
%fluid 1 boundary
[ F1 V1 ] = isosurface(g.xs{1},g.xs{2},g.xs{3},data1,isolevel);
v1_sz = size(V1); v1_sz = v1_sz(1)
f1_sz = size(F1); f1_sz = f1_sz(1);

% fluid 2 boundary
[ F2 V2 ] = isosurface(g.xs{1},g.xs{2},g.xs{3},data2,isolevel);
v2_sz = size(V2); v2_sz = v2_sz(1)
f2_sz = size(F2); f2_sz = f2_sz(1);

%find vertices and triangles on fluid1-fluid2 interface
[ V_FLD_INTFC i1 i2] = intersect(V1,V2,'rows');

%find triangles that have all 3 vertices in V_FLD
v1_fld = zeros(v1_sz,1); v1_fld(i1) = 1;
f_fld = find((v1_fld(F1(:,1)) & v1_fld(F1(:,2)) & v1_fld(F1(:,3))));
F_FLD_INTFC = F1(f_fld,:);

%isolate all traingles that are not in F_FLD - for visualization purpose
f_grain = find( (v1_fld(F1(:,1))==0) | (v1_fld(F1(:,2))==0) | (v1_fld(F1(:,3))==0) );
F_GRAIN = F1(f_grain,:);

% visualize
if( do_vis )
    color = 'y'; %yellow
    
    figure, h_fld = patch('Faces',F_FLD_INTFC,'Vertices',V1);
    set(h_fld,'FaceColor',color,'EdgeColor', 'none');
    daspect([1 1 1]); camlight; lighting phong;
    view(3), hold on
   
    h_grain = patch('Faces',F_GRAIN,'Vertices',V1);
    size(F_GRAIN)
    set(h_grain, 'FaceColor', [0.8 0.8 0.8] ,'EdgeColor', 'none'); % gray color
    daspect([1 1 1]); camlight; lighting phong;
    view(3), hold off
end

%find triple points
V = find_bdry_vertices(V1,F_FLD_INTFC);

dx_mask = centeredFirstSecond(g, mask, 1);
dy_mask = centeredFirstSecond(g, mask, 2);
dz_mask = centeredFirstSecond(g, mask, 3);

dx_data1= centeredFirstSecond(g, data1, 1);
dy_data1= centeredFirstSecond(g, data1, 2);
dz_data1= centeredFirstSecond(g, data1, 3);

%interpolate values at points of interest
dx_mask = interp3(dx_mask,V(:,1),V(:,2),V(:,3));
dy_mask = interp3(dy_mask,V(:,1),V(:,2),V(:,3));
dz_mask = interp3(dz_mask,V(:,1),V(:,2),V(:,3));
 
dx_data1= interp3(dx_data1,V(:,1),V(:,2),V(:,3));
dy_data1= interp3(dy_data1,V(:,1),V(:,2),V(:,3));
dz_data1= interp3(dz_data1,V(:,1),V(:,2),V(:,3));

prod = dx_mask .* dx_mask + dy_mask .* dy_mask + dz_mask .* dz_mask;
prod1 = dx_data1 .* dx_data1 + dy_data1 .* dy_data1 + dz_data1 .* dz_data1;

prod = prod*prod1;  prod = sqrt(prod);
prod( find( prod == 0) ) = 1; % avoid division by zero

prod3 = dx_mask .* dx_data1 + dy_mask .* dy_data1 + dz_mask .* dz_data1;
cos_angle = prod3 ./ prod; 
% results > 1.0 or < -1.0 are rounding off errors, correct for those
cos_angle( find( cos_angle > 1.0) ) = 1.0;
cos_angle( find( cos_angle < -1.0) ) = -1.0;

%evaluate cos_angle on triple point vertices
cos_angle_intrp = interp3(cos_angle,V(:,1),V(:,2),V(:,3));
angle_intrp = acos(cos_angle_intrp);
angle_intrp = angle_intrp * (180 / pi);

if( do_vis ) figure, hist(cos_angle,50); end;

num_of_triple_points = size(cos_angle)
median_angle  = median(cos_angle)


% normalize gradient vectors to get outward normal for plot
dx_mask = dx_mask ./ prod; dy_mask = dy_mask ./ prod;   dz_mask = dz_mask ./ prod;
dx_data1=dx_data1 ./ prod1; dy_data1=dy_data1 ./ prod1; dz_data1=dz_data1 ./ prod1;

[X Y Z] = meshgrid(V(:,1),V(:,2),V(:,3));
size(X)
dxm = interp3(dx_mask,V(:,1),V(:,2),V(:,3));
dym = interp3(dy_mask,V(:,1),V(:,2),V(:,3));
dzm = interp3(dz_mask,V(:,1),V(:,2),V(:,3));
[U V W] = meshgrid(dxm,dym,dzm);

figure, quiver3(X,Y,Z,U,V,W);

if( do_vis )
    figure, h_fld = patch('Faces',F_FLD_INTFC,'Vertices',V1);
    set(h_fld,'FaceColor',color,'EdgeColor', 'none');
    daspect([1 1 1]); camlight; lighting phong;
    view(3), hold on
    
    disp('Values of contact angle below median are plot in blue, otherwise magenta.');
    cut_off1 = median_angle; cut_off2 = median_angle;
    i = find( angle_intrp <= cut_off1);
    j = find( (angle_intrp > cut_off1 )  & (angle_intrp < cut_off2 ) );
    k = find(angle_intrp > cut_off2);

    plot3(V(i,1),V(i,2),V(i,3),'b.')
    plot3(V(j,1),V(j,2),V(j,3),'c.')
    plot3(V(k,1),V(k,2),V(k,3),'m.')
end
