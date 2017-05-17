function opposite_side = findOppSide(mask,g,data0)

%  function opposite_side = findOppSide(mask,g,data0)
%  opposite_side are voxel indices that can be checked later on if reached
%  by an invading fluid (if data(opposite_side) < 0 for at least one
%  voxel).
%  These voxels are pore voxels with exclusion of voxels occupied by fluid
%  at initial time (data0).

opposite_side = [];
if(g.dim ~= 2) 
    disp('findOppSide() works in 2D only.')
    return;
end

side_xl = find(mask(1,:) < 0); % pore voxels on lower x side
ic_xl =  find(data0(1,:) < 0); % initial position fluid voxels on lower x side
side_xl = setdiff(side_xl,ic_xl);

if( numel(side_xl) ) 
    i = ones(size(side_xl));
    opposite_side = sub2ind(size(mask),i,side_xl);
end


side_xu = find(mask(g.N(1),:) < 0);
ic_xu = find(data0(g.N(1),:) < 0);
side_xu = setdiff(side_xu,ic_xu);

if( numel(side_xu) ) 
    i = ones(size(side_xu))*g.N(1);
    opposite_side = union(opposite_side,sub2ind(size(mask),i,side_xu));
end

side_yl = find(mask(:,1) < 0);
ic_yl = find(data0(:,1) < 0);
side_yl = setdiff(side_yl,ic_yl);

if( numel(side_yl) ) 
    j = ones(size(side_yl));
    opposite_side = union(opposite_side,sub2ind(size(mask),side_yl,j));
end

side_yu = find(mask(:,g.N(2)) < 0);
ic_yu = find(data0(:,g.N(2)) < 0);
side_yu = setdiff(side_yu,ic_yu);

if( numel(side_yu) ) 
    j = ones(size(side_yu))*g.N(2);
    opposite_side = union(opposite_side,sub2ind(size(mask),side_yu,j));
end