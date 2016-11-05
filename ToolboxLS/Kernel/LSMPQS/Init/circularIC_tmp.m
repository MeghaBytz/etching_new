function data = circularIC_tmp(g,mask)
% data = circularICx(g)
% create circular front as the initial condition
% g - grid

maskCenter = [g.min(1); (g.min(2) + g.max(2))*0.5;0; 0.0; 0.0 ];
tmp = size(find(mask(1,:) < 0));
tmp = tmp(1) *tmp(2); %initial slab of voxels on entry
maskRadius = (tmp/2.0 + 1)*g.dx(1);

data =  shapeSphere(g, maskCenter, maskRadius);

%this is for half-circle
normal = [ 1.0; 0.0];
point  = [(g.min(1) + g.max(1))*0.5;g.min(2)];
data = min(data, shapeHyperplane(g,normal,point));
    