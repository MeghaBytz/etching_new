function data = circularICx(g,mask)
% data = circularICx(g)
% create circular front as the initial condition
% g - grid

pos = 1;
[i j] = find(mask(pos,:) < 0); %coords initial slab of voxels on entry
tmp = numel(j);
mid = floor(tmp/2);
maskCenter = [g.min(1); g.min(2) + j(mid)*g.dx(2);0; 0.0; 0.0 ]; %midpoint of initial slab
maskRadius = (tmp/2.0 + 1)*g.dx(1);

data =  shapeSphere(g, maskCenter, maskRadius);

% this is for half-circle
%normal = [ 1.0; 0.0];
%point  = [(g.min(1) + g.max(1))*0.5;g.min(2)];
%data = min(data, shapeHyperplane(g,normal,point));
    