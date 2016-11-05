function [curv_centers, heights] = getBottomCurvatures(g, data, layer_centers, init_level)

center_ind = zeros(length(layer_centers),1);

% get layer center x-indices
for i = 1:length(layer_centers)
    center_ind(i) = find(abs(g.xs{1,1}(:,1) - layer_centers(i)) < g.dx(1));
end

% calculate curvature everywhere
curv = curvatureSecond(g,data);

% find interface points
intfcPoints = isNearInterface(data);

% initialize heights
h = zeros(size(data,2),1);
% get bottom y-indices
bottom_y_ind = zeros(size(data,2),1);

% locate index of initial level
% top_ind = find(abs(g.xs{2,1}(1,:)-init_level) < g.dx(1));
top_ind = 1;
% find heights at each x
for i = 1:size(data,1)
    %top_end = 0; bottom_end = 0;
    %bottom_ind = 0;
    ind = find(intfcPoints(i,:));
    
    %if(size(ind,2) == 2)
        bottom_ind = ind(2);
        bottom_end = data(i,ind(2));
        top_end = 0;
    %end
    %h(i) = (top_ind  - bottom_ind)*g.dx(1) + top_end - bottom_end;
    h(i) = (bottom_ind - top_ind)*g.dx(1) - bottom_end;
    bottom_y_ind(i) = bottom_ind;
end


% get heights only at layer centers
heights = h(center_ind);

% get bottom ends of the trenches at layer centers
bottoms = bottom_y_ind(center_ind);

%remove first and last layers
center_ind(1) = []; center_ind(end) = [];
bottoms(1) = []; bottoms(end) = [];
heights(1) = []; heights(end) = [];

curv_centers = zeros(length(center_ind),1);
% get curvatures knowing the layer centers and the bottom of the trenches
for i = 1:length(center_ind)
    curv_centers(i) = curv(sub2ind(size(data),center_ind(i), bottoms(i)));
end

end