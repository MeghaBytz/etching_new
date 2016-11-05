function [startMatInd, endMatInd, map] = getMaterialMap(g,options)
%define start x and y indices for each material based on -1 to 1 grid
%test case shown here for a row of materials
%Algernon/Kepler/Algernon/Kepler/Algernon
%{
startXIndices = [g.min(1)];
endXIndices = [g.max(1)];
startYIndices = [g.min(2)];
endYIndices = [g.max(2)];
startMatInd = [startXIndices;startYIndices]
endMatInd = [endXIndices;endYIndices]
centMatInd = [(startXIndices + endXIndices)/2];
map = 1;
%}

if(options.horizontal) % horizontal stacks
    startXIndices = g.min(1)*[1 1 1 1 1];
    endXIndices = g.max(1)*[1 1 1 1 1];
    startYIndices = [-1 -.6 -.2 .2 .6];
    endYIndices = [-.6 -.2  .2 .6 1];
    startMatInd = [startXIndices;startYIndices];
    endMatInd = [endXIndices;endYIndices];
    map = [1 2 1 2 1];    
else                    % vertical stacks
    startXIndices = [-1 -.6 -.2 .2 .6];
    endXIndices = [-.6 -.2  .2 .6 1];
    startYIndices = [1 1 1 1 1];
    endYIndices = [-1 -1 -1 -1 -1];
    startMatInd = [startXIndices;startYIndices];
    endMatInd = [endXIndices;endYIndices];
    map = [1 2 1 2 1];
end
