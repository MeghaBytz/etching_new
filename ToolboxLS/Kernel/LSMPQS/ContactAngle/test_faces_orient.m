
function [orient NF] = test_faces_orient(F,V,data)

% Test triangle normals' orientation
% F,V, data  - faces, vertices and data they were calculated from using
% isosurface
% orient(f,i) is a dot product of face f normal and OUTWARD vertex normal at it's
% vertex i, i=1,2,3
% NF(f,:) is the triangle f normal

NV = isonormals(data,V,'negate'); %point OUTWARD, to the larger 'data' values
v_sz = size(V); v_sz = v_sz(1);
f_sz = size(F); f_sz = f_sz(1);

for(f=1:f_sz)
    v(1) = F(f,1);  v(2) = F(f,2); v(3) = F(f,3);
    p = V(v,:);
    % triangle normal, not normalized
    n = cross(p(2,:)-p(1,:),p(3,:)-p(1,:));
    nnorm = norm(n); if (nnorm == 0) nnorm = 1; end
    NF(f,:) = n/nnorm;
    
    % dot product of triangle normal with the normal at each vertex
    orient(f,1) = dot(n,NV(v(1),:));
    orient(f,2) = dot(n,NV(v(2),:));
    orient(f,3) = dot(n,NV(v(3),:));
end
