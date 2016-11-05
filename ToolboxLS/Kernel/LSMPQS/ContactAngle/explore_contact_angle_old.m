function  [V angle_intrp] = explore_contact_angle(infl_name1,infl_name2)

% [V angle_intrp] = explore_contact_angle(infl_name1,infl_name2)
% explores contact angle btw fluid1-fluid2 interface and fluid1-grain
% interface
% infl_name1 - filename for file containing the signed distance function of
%               fluid blob of interest (produced by 3DMA-Rock)
% infl_name2 - filename with signed distance function for the complimentary
%               fluid phase in the blob subvolume (produced by 3DMA-Rock)
% V - vertices from fluid1-fluid2 interface surface identified as triple
% points
% angle_interp = contact angle interpolated at V

run('~/ToolboxLS-1.0/Examples/addPathToKernel.m');
addpath ~/3dma_rock/src/matlab  %for read_l2_brnfl

% data1 < 0 defines region1 occupied by fluid1 and similarly data2
[data1 nx ny nz] = read_l2_brnfl(infl_name1);
[data2 nx ny nz] = read_l2_brnfl(infl_name2);

mask = min(data1,data2);  %pore space level set function, commonly referred to as mask

% create grid, needed for curvature computation (with LS Toolbox function)
g.dim = 3;
g.min = 1;
g.dx = 1; % irrelevant for now
g.max = g.dx*[nx; ny; nz];
g.bdry = @addGhostExtrapolate;
g = processGrid(g)

isolevel = 0;

% fluid 1 boundary
[ F1 V1 ] = isosurface(g.xs{1},g.xs{2},g.xs{3},data1,isolevel);
v1_sz = size(V1); v1_sz = v1_sz(1)
f1_sz = size(F1); f1_sz = f1_sz(1);

% fluid 2 boundary
[ F2 V2 ] = isosurface(g.xs{1},g.xs{2},g.xs{3},data2,isolevel);
v2_sz = size(V2); v2_sz = v2_sz(1)
f2_sz = size(F2); f2_sz = f2_sz(1);

%find vertices on fluid1-fluid2 interface
FLD_INTFC = intersect(V1,V2,'rows');

%find triple points (V are neighboring fluid1 & fluid2, looking for ones
%that neighbor grain as well
V = findEdgeVertex3D(FLD_INTFC,mask); 

dx_mask = centeredFirstSecond(g, mask, 1);
dy_mask = centeredFirstSecond(g, mask, 2);
dz_mask = centeredFirstSecond(g, mask, 3);

dx_data1= centeredFirstSecond(g, data1, 1);
dy_data1= centeredFirstSecond(g, data1, 2);
dz_data1= centeredFirstSecond(g, data1, 3);

prod = dx_mask .* dx_mask + dy_mask .* dy_mask + dz_mask .* dz_mask;
prod = sqrt(prod);

prod1 = dx_data1 .* dx_data1 + dy_data1 .* dy_data1 + dz_data1 .* dz_data1;
prod1 = sqrt(prod1);

prod3 = dx_mask .* dx_data1 + dy_mask .* dy_data1 + dz_mask .* dz_data1;

% avoid division by zero when one of vectors is a zero vector
prod( find( prod == 0) ) = 1;
prod1( find( prod1 == 0) ) = 1;

cos_angle = prod3 ./ (prod .* prod1); 
find( isnan(cos_angle) == 1 )
% results > 1.0 or < -1.0 are rounding off errors, correct for those
cos_angle( find( cos_angle > 1.0) ) = 1.0;
cos_angle( find( cos_angle < -1.0) ) = -1.0;

%evaluate cos_angle on triple point vertices
extrapol_val = 1.0 % set extrapolated points to 1.0
cos_angle_intrp = interp3(cos_angle,V(:,1),V(:,2),V(:,3), 'linear',extrapol_val);
angle_intrp = acos(cos_angle_intrp);
angle_intrp = angle_intrp * (180 / pi);
hist(angle_intrp,50);
size(angle_intrp)
median_angle  = median(angle_intrp)
size(find( abs(angle_intrp) < 5))


function  Vlist = findEdgeVertex3D(FLD_INTFC,mask)

%  findEdgePoint3D(FLD_INTFC,mask)

%  Finds fluid interface vertices that have grain phase voxels in their
%  neighborhood.
%  FLD_INTFC = vertices from triangulated surface between two fluids
%              note that each row specifies a vertex
%  mask < 0 defines where pore space is (mask > 0 defines grain space)
%  g -  grid

fld_sz = size(FLD_INTFC); fld_sz = fld_sz(1);

% for each vertex (row inFLD_INTFC) only one index is on a edge between
% voxels from different fluids. e.g. v = (4.5,6,7) identifies vertex
% between voxels A(4,6,7) and B(5,6,7). We need to identify points A and B
% and test their neighbors

FLD_floor = floor(FLD_INTFC); % basically all 'A' points
FLD_ceil = ceil(FLD_INTFC);  % basically all 'B' points

[nx ny nz] = size(mask);
 
count = 0; % counter of edge vertices

for(k=1:fld_sz)
   test = 0; % test if grain voxel found in nbhd
   
   A = [ FLD_floor(k,1) FLD_floor(k,2) FLD_floor(k,3) ];
   A = max([1 1 1], A);
   A = min([nx ny nz],A);
   
   B = [ FLD_ceil(k,1)  FLD_ceil(k,2) FLD_ceil(k,3) ];
   B = max([1 1 1], B);
   B = min([nx ny nz],B);
   
   test =  test_nbhd(A,mask,nx,ny,nz);
   if( test )
       count = count + 1;
       Vlist(count,:) = FLD_INTFC(k,:);
   else
       test =  test_nbhd(B,mask,nx,ny,nz);
       if( test )
          count = count + 1;
          Vlist(count,:) = FLD_INTFC(k,:);
       end
   end
end

count

function test =  test_nbhd(A,mask,nx,ny,nz)

% tests 6-neighborhood of voxel A for grain phase
test = 0;

xu = min(A(1)+1,nx);
if( mask(xu,A(2),A(3)) > -10*eps ) 
    test = 1; return;
end

xl = max(A(1)-1,1);
if( mask(xl,A(2),A(3)) > -10*eps ) 
    test = 1; return;
end

yu = min(A(2)+1,ny);
if( mask(A(1),yu,A(3)) > -10*eps ) 
    test = 1; return;
end

yl = max(A(2)-1,1);
if( mask(A(1),yl,A(3)) > -10*eps ) 
    test = 1; return;
end

zu = min(A(3)+1,nz);
if( mask(A(1),A(2),zu) > -10*eps ) 
    test = 1; return;
end

zl = max(A(3)-1,1);
if( mask(A(1),A(2),zl) > -10*eps ) 
    test = 1; return;
end
