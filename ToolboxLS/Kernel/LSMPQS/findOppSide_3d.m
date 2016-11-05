function opposite_side = findOppSide_3d(mask,g,data0)

%  function opposite_side = findOppSide(mask,g,data0)
%  opposite_side are voxel indices that can be checked later on if reached
%  by an invading fluid (if data(opposite_side) < 0 for at least one
%  voxel).
%  These voxels are pore voxels with exclusion of voxels occupied by fluid
%  at initial time (data0).

opposite_side = [];

side_zl = find(mask(:,:,1) < 0); % pore voxels on lower z side (entry side)
ic_zl =  find(data0(:,:,1) < 0); % initial position fluid voxels on lower z side
side_zl = setdiff(side_zl,ic_zl);

if( numel(side_zl) ) 
    i = ones(size(side_zl));
    opposite_side = sub2ind(size(mask),i,side_zl);
end


side_zu = find(mask(:,:,g.N(3)) < 0);
ic_zu = find(data0(:,:,g.N(3)) < 0);
side_zu = setdiff(side_zu,ic_zu);

if( numel(side_zu) ) 
    i = ones(size(side_zu))*g.N(3);
    opposite_side = union(opposite_side,sub2ind(size(mask),i,side_zu));
end

%{

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
%}
end
