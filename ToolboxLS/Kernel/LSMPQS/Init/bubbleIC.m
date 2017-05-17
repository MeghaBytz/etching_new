function data = bubbleIC(g)

% data = bubbleIC(g)
% create circular front as the initial condition
% g - grid

maskCenter = [(g.min(1) + g.max(1))*0.5; 0.0; 0.0; 0.0 ];
maskRadius = 0.5;

data =  shapeSphere(g, maskCenter, maskRadius);

% this is for half-circle
%normal = [ 1.0; 0.0];
%point  = [(g.min(1) + g.max(1))*0.5;g.min(2)];
%data = min(data, shapeHyperplane(g,normal,point));
    