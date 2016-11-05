function plotLevelSetInterior(data, level, mask)

% function plotLevelSetInterior(data, level, mask)
% Function plots points that satisfy 'data < level' in white, 'mask >
% level' in black (masked region).

if( nargin < 2 ) 
    error('plotLevelSetInterior needs at least 2 arguments');
end

plotData = 2*ones(size(data));
if( nargin == 3) 
    plotData( find(mask > level) ) = 0;
end

plotData( find(data < level) ) = 1;
pcolor(plotData'); colormap('hot'); axis equal; axis tight;