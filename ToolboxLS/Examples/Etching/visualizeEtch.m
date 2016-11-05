function h = visualizeEtch(g, data, mask, options, title_str)

if(options.horizontal)  % layers are arranged horizontally
    layers = [options.layer_boundaries(2,:); options.layer_boundaries(4,:)]';
else                    % layers are arranged vertically
    layers = [options.layer_boundaries(1,:); options.layer_boundaries(3,:)]';
end

colors = options.colors;
nLayers = size(layers,1);
layer_colors = zeros(nLayers, 3);
data_plot_1 = data; data_plot_2 = data; data_plot_3 = data;
 
for i = 1:length(colors)
    switch colors(i)
        case 'b'
            layer_colors(i,:) = [0.1905    0.8095    1.0000];
        case 'r'
            layer_colors(i,:) = [0.6875 0 0];
        case 'y'
            layer_colors(i,:) = [0.9955 0.7861 0.1967];
        case 'g'
            layer_colors(i,:) = [];
        case 'p'
            layer_colors(i,:) = [];
    end
end

% Apply colors for each layer
for i = 1:nLayers
    if(options.horizontal)
        data_plot_1((g.xs{2,1} >= layers(i,1)) & g.xs{2,1} <= layers(i,2)) = layer_colors(i,1);
        data_plot_2((g.xs{2,1} >= layers(i,1)) & g.xs{2,1} <= layers(i,2)) = layer_colors(i,2);
        data_plot_3((g.xs{2,1} >= layers(i,1)) & g.xs{2,1} <= layers(i,2)) = layer_colors(i,3);
    else
        data_plot_1((g.xs{1,1} >= layers(i,1)) & g.xs{1,1} <= layers(i,2)) = layer_colors(i,1);
        data_plot_2((g.xs{1,1} >= layers(i,1)) & g.xs{1,1} <= layers(i,2)) = layer_colors(i,2);
        data_plot_3((g.xs{1,1} >= layers(i,1)) & g.xs{1,1} <= layers(i,2)) = layer_colors(i,3);
    end
end
 
% White (blank) above everything
data_plot_1(data >= 0) = 1; data_plot_2(data >= 0) = 1; data_plot_3(data >= 0) = 1;

% Apply mask (black)
if(options.doMask)
    data_plot_1(mask >= 0) = 0; data_plot_2(mask >= 0) = 0; data_plot_3(mask >= 0) = 0;
end

%concatenate everything into a single colored image and image it
data_plot = cat(3, data_plot_1, data_plot_2, data_plot_3);

% crop out the boundaries
% find number of pixels to crop out - should be half the number of x-points
% in the boundary layers - should be necessary only in vertical layers case
if(~options.horizontal)
    cropRight = (g.xs{1,1} <= layers(1,1) & g.xs{1,1} >= 0.5*(layers(1,2) + layers(1,1)));
    nPointsRight = size(find(cropRight(:,1)),1);
    cropLeft = (g.xs{1,1} <= layers(end,1) & g.xs{1,1} >= 0.5*(layers(end,2) + layers(end,1)));
    nPointsLeft = size(find(cropLeft(:,1)),1);
    crop_data_plot = imcrop(rot90(data_plot),[nPointsLeft 0 (size(data,1)-nPointsRight-nPointsLeft) size(data,2)]);
    h = image(crop_data_plot);
else
    h = image(rot90(data_plot));
end

title(title_str);

assignin('base', 'dataFromLevelSEt', data)
assignin('base', 'dataMapped', data_plot)
end