function [a,b,vel] = abplot_new_1(mask,delx,a0,ca, g)

% function [a,b,vel] = abplot(mask,delx,a0,ca)
% 1. This function generates matrices to input into the level set equation:
% phi_t +(a-b*kappa)|del phi| + v.(del phi) = 0
% The a, b and velocity terms vary spatially and hence need to be
% implemented as matrices. The function takes as input:
% mask: geometry of the domain
% delx: grid size, equal for all directions
% a0: pressure like term
% ca: contact angle imposed

% 4. test accuracy - medium vs higher, also the value of epsilon for
% heaviside function

% 7. can handle 3d grids.
% 
gradx = centeredFirstSecond(g, mask, 1);
grady = centeredFirstSecond(g, mask, 2);
if(g.dim>2)
    gradz = centeredFirstSecond(g,mask,3);
end

if(g.dim<3)
    absgrad = (gradx.^2+grady.^2).^0.5;
else
    absgrad = (gradx.^2+grady.^2+gradz.^2).^0.5;
end

C1 = 0.04;
C2 = 1;%1/delx;
h1 = h_side(-1*mask,delx);
h2 = h_side(mask,delx);

b = h1*C1;
k0 = a0;

s = signfunc(mask, absgrad, delx);

vel{1,1} = s.*h2.*gradx*C2;
vel{2,1} = s.*h2.*grady*C2;
if(g.dim>2)
    vel{3,1} = s.*h2.*gradz*C2;
end
a = h1.*k0-s.*h2.*absgrad.*cos(ca)*C2;

