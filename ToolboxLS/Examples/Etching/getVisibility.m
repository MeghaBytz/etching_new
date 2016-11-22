function vis = getVisibility(g, source, data, mask)

% get points where data > 0
pts = find((data>0) && (mask < 0));
d = zeros(size(pts));
vis = zeros(size(data));
% iterate through pts
for i = pts
    cur_ind = ind2sub(size(data),i);
    % draw line segment through pt i and source
    a = source - cur_ind;
    % at each point in domain, find normal distance
    % if normal distance is less than grid dx, and level set value is less
    % than 0, points can't see each other
    for j = 1:size(data(:))
        b_ind = ind2sub(size(data),j);
        b = b_ind - cur_ind;
        d(j) = norm(cross(a,b))/norm(a);
        if (d(j) <= g.dx(1) && data(j)<0)
            vis(i) = 0;
        end 
    end
end
