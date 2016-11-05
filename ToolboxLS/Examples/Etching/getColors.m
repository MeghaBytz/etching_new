function colors = getColors(map)

% define available colors for layers
available_colors = ['b'; 'r'; 'y'; 'g'; 'p'];

colors = num2str(zeros(size(map,2),1));

% assign colors to each layer according to the material map
for i = 1:size(map,2)
    colors(i) = available_colors(map(i));
end

end