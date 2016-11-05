function  [mask,g] = maskCorner(g,p)

% function  [mask,g] = maskCorner(g)
% Create mask that corresponds to a corner formed by three
% points p(1), p(2), p(3) or rather lines (p(1),p(2)) and (p(2),p(3))

mask = zeros(g.shape);

for(i=1:2)
    %find normal of the line p(i+1),p(i)
    v  = p(:,i+1) - p(:,i);
    normv = norm( v );
    v = v/normv;
    n = [v(2); -v(1) ];

    mask1 = shapeHyperplane(g,-n,p(:,i));
    mask = min(mask,mask1);
end
