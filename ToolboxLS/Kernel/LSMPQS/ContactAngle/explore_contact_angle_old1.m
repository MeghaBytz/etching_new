function  [median_angle angle_intrp] = explore_contact_angle(infl_name1,infl_name2,do_vis)

% [V angle_intrp] = explore_contact_angle(infl_name1,infl_name2,do_vis)
% explores contact angle btw fluid1-fluid2 interface and fluid1-grain
% interface
% infl_name1 - filename for file containing the signed distance function of
%               fluid blob of interest (produced by 3DMA-Rock)
% infl_name2 - filename with signed distance function for the complimentary
%               fluid phase in the blob subvolume (produced by 3DMA-Rock)
% do_vis - set to 1 if you want to visualaize result
% V - vertices from fluid1-fluid2 interface surface identified as triple
% points
% angle_interp = contact angle interpolated at V

if (nargin < 3) do_vis = 0; end
    
run('~/ToolboxLS-1.0/Examples/addPathToKernel.m');
addpath ~/3dma_rock/src/matlab  %for read_l2_brnfl

% data1 < 0 defines region1 occupied by fluid1 and similarly data2
[data1 nx ny nz] = read_l2_brnfl(infl_name1);
[data2 nx ny nz] = read_l2_brnfl(infl_name2);

mask = min(data1,data2);  %pore space level set function
mask = -mask;             %grain space level set function

% create grid, needed for curvature computation (with LS Toolbox function)
g.dim = 3;
g.min = 1;
g.dx = 1; % irrelevant for now
g.max = g.dx*[nx; ny; nz];
g.bdry = @addGhostExtrapolate;
g = processGrid(g)

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
    
    if( findstr('H2O',infl_name1) )     color = [0 0.8 0]; %green
    elseif( findstr('OIL',infl_name1) ) color = [0.8 0 0]; %red
    else                                color = 'y'; %yellow
    end
    
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

%find vertices on the fluid1-fluid2 interface boundary
V = find_bdry_vertices(V1,F_FLD_INTFC);

dx_mask = centeredFirstSecond(g, mask, 1);
dy_mask = centeredFirstSecond(g, mask, 2);
dz_mask = centeredFirstSecond(g, mask, 3);

dx_data1= centeredFirstSecond(g, data1, 1);
dy_data1= centeredFirstSecond(g, data1, 2);
dz_data1= centeredFirstSecond(g, data1, 3);

dx_data2= centeredFirstSecond(g, data2, 1);
dy_data2= centeredFirstSecond(g, data2, 2);
dz_data2= centeredFirstSecond(g, data2, 3);

%interpolate values at points of interest
method = 'linear';
dxm = interp3(dx_mask,V(:,1),V(:,2),V(:,3),method);
dym = interp3(dy_mask,V(:,1),V(:,2),V(:,3),method);
dzm = interp3(dz_mask,V(:,1),V(:,2),V(:,3),method);
 
dxd1= interp3(dx_data1,V(:,1),V(:,2),V(:,3),method);
dyd1= interp3(dy_data1,V(:,1),V(:,2),V(:,3),method);
dzd1= interp3(dz_data1,V(:,1),V(:,2),V(:,3),method);

dxd2= interp3(dx_data2,V(:,1),V(:,2),V(:,3),method);
dyd2= interp3(dy_data2,V(:,1),V(:,2),V(:,3),method);
dzd2= interp3(dz_data2,V(:,1),V(:,2),V(:,3),method);

fprintf('Calculation through the first fluid, %s',infl_name1);
angle1 = calculate_angle(dxm,dym,dzm,dxd1,dyd1,dzd1,do_vis);
fprintf('Calculation through the second fluid, %s',infl_name2);
angle2 = calculate_angle(dxm,dym,dzm,dxd2,dyd2,dzd2,do_vis);

fprintf('dx_data1 and dx_data2 angle on the boundary');
calculate_angle(dxd1,dyd1,dzd1,dxd2,dyd2,dzd2,do_vis);

plot_velocity_Geomview('tmp_mask.list',V(:,1),V(:,2),V(:,3),dx_mask,dy_mask,dz_mask,[0.0 0.0 0.0],g.min,g.max);
plot_velocity_Geomview('tmp_data1.list',V(:,1),V(:,2),V(:,3),dxd1,dyd1,dzd1,[0.0 0.0 0.5],g.min,g.max);
plot_velocity_Geomview('tmp_data2.list',V(:,1),V(:,2),V(:,3),dxd2,dyd2,dzd2,[0.5 0.5 0.0],g.min,g.max);

%interpolate at all fld1-fld2 interface points
method = 'linear';
Vint = setdiff(V_FLD_INTFC,V,'rows');

%x = V_FLD_INTFC(:,1); y = V_FLD_INTFC(:,2);   z = V_FLD_INTFC(:,3);
x = Vint(:,1);  y = Vint(:,2);   z = Vint(:,3);
dxd1= interp3(dx_data1,x,y,z,method);
dyd1= interp3(dy_data1,x,y,z,method);
dzd1= interp3(dz_data1,x,y,z,method);

dxd2= interp3(dx_data2,x,y,z,method);
dyd2= interp3(dy_data2,x,y,z,method);
dzd2= interp3(dz_data2,x,y,z,method);

fprintf('dx_data1 and dx_data2 angle on interior fluid interface');
calculate_angle(dxd1,dyd1,dzd1,dxd2,dyd2,dzd2,do_vis);

plot_velocity_Geomview('tmp_data1_all.list',x,y,z,dxd1,dyd1,dzd1,[0.0 0.0 0.5],g.min,g.max);
plot_velocity_Geomview('tmp_data2_all.list',x,y,z,dxd2,dyd2,dzd2,[0.5 0.5 0.0],g.min,g.max);

if( do_vis )
    figure, h_fld = patch('Faces',F_FLD_INTFC,'Vertices',V1);
    set(h_fld,'FaceColor',color,'EdgeColor', 'none');
    daspect([1 1 1]); camlight; lighting phong;
    view(3), hold on
    
    disp('Values of contact angle below median are plot in blue, otherwise magenta.');
    median_angle = median(angle1);
    cut_off1 = median_angle; cut_off2 = median_angle;
    i = find( angle1 < cut_off1);
    j = find( (angle1> cut_off1 )  & (angle1 < cut_off2 ) );
    k = find(angle1> cut_off2);

    plot3(V(i,1),V(i,2),V(i,3),'b.')
    plot3(V(j,1),V(j,2),V(j,3),'c.')
    plot3(V(k,1),V(k,2),V(k,3),'m.')
end