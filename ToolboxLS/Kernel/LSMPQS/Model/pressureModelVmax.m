function new_a = pressureModelVmax(t,data,schemeData)

%--------------------------------------------------------------------------
%  Model for varying velocity in normal direction
%--------------------------------------------------------------------------
% function new_a = pressureModelVmax(t,data,schemeData)
%
% a(t)=new_a is capillary pressure modeled as 
% a(t) = a0 * exp(-f*(V(t)- Vmax)/Vmax) where V(t) is 
% volume of the invading fluid at time t, and Vmax the maximal volume the
% fluid could occupy.
% a(t) --> new_a
% f --> schemeData.factor
% a0 --> schemeData.magnitude
% Vmax --> schemeData.vmax
% V(t) volume of the set described as data<0

if( schemeData.grid.dim == 3 )
    volume = size(find(data < 0));
    volume = volume(1); %find current volume occupied by invading fluid
else
   volume = areaLevelSetInterior(data,schemeData.grid);
end

% change relative to the max. possible volume
vmax = schemeData.vmax; % max possible volume - kind-a arbitrary
rel_change = (volume - vmax)/vmax;
factor = schemeData.factor;  % arbitrary, set beforehand
new_a = schemeData.magnitude * exp(-factor*rel_change);
%schemeData.v0 = volume;