function area = areaLevelSetInterior(data,g)

% function area = areaLevelSetInterior(data,g)
% data - level set function
% g - grid
% area - (polygonal) area for the region data < 0 (interior of level set func)
% This measure should be much more accurate than simply finding number of
% points such that data < 0

 sz = size( find(data < 0)); sz = sz(1); % rough area measure

 if ( sz == 0 )
     area = 0;
     return;
 end
 
 % find a 0-level set and output coords to matrix C 
 C = contourc(data,[0 0]);
 % get the coordinates of all contour points
 [C1 n]= getContourPoints(C);
 x = C1(1,:);
 y = C1(2,:);
 
 area = polyarea(y,x);

 % if area is more than 20% different than 'sz' there might be an issue
 % e.g. if contour encloses disjoint sets, polyarea is no good
 r = abs( (area - sz)/sz );
 if ( abs( (area - sz)/sz ) > 0.2)
     %fprintf('\nsz %d area %g  - ',sz,area);
     %fprintf('Polygonal area replaced by pixel count in areaLevelSetInterior()');
     area = sz;
 end