
function [distance,x11, x12] = anSolCirclesGeneral(R,theta,R1,x,y,y1)

% function [distance,x11, x12] = anSolCirclesGeneral(R,theta,R1,x,y,y1)
% R,R1 - circle radii
% theta - angle in degrees
%
% returns distance between centers of two circles of radii R and R1
% such that they cross at angle theta
% if (x,y) coordinate of the first center are given, and y1 of the second
% then x1 is calculated (there are two possibilities)

% We have to solve a triangle (a,b,c,alpha,beta,gamma) where we have given
% two sides a,b=R,R1 that enclose gamma=180 - theta. The third side is the
% distance. Draw to see.

gamma = 180 - theta
a = R
b = R1

% we know alpha+beta = 180-gamma=theta
tmp = cotd(0.5*gamma)*(a-b)/(a+b); % tan((alpha-beta)/2)
tmp1 = 2*atand(tmp); %alpha-beta
tmp2 = theta;        %alpha+beta
alpha = 0.5* ( tmp1 + tmp2)
beta = theta-alpha

c=a*sind(gamma)/sind(alpha);

if( (alpha >= 0) && (beta >= 0)) 
    distance=c;
else 
     fprintf('Error: alpha %g beta %g',alpha,beta);
     distance = -1;
end

if (nargin >= 6)
    xdiff_sq = distance*distance - (y-y1)*(y-y1);
    x11 = x + sqrt(xdiff_sq);
    x12 = x - sqrt(xdiff_sq);
end    
return;



