function [a,b,vel] = abplot_new(mask,delx,a0,ca, g)

% function [a,b,vel] = abplot(mask,delx,a0,ca)
% 1. COMMENT WHAT INPUTS AND OUTPUTS ARE
%This function generates matrices to input into the level set equation:
% phi_t +(a-b*kappa)|del phi| + v.(del phi) = 0
% The a, b and velocity terms vary spatially and hence need to be
% implemented as matrices. The function takes as input:
% mask: geometry of the domain
% delx: grid size, equal for all directions
% a0: pressure like term
% ca: contact angle imposed


% 2. have an accompanying file mentioning the equation modification for our
% code
% 3. testing -- make sure what is the range of angles where it works well
% and how does it compare to the WRR paper
% 
% 4. test accuracy - medium vs higher, also the value of epsilon for
% heaviside function

% 5. avoid for loops
% 6. Check reinitialization
% 
gradx = centeredFirstSecond(g, mask, 1);
grady = centeredFirstSecond(g, mask, 2);
absgrad = (gradx.^2+grady.^2).^0.5;
[m, n] = size(mask);

C1 = 0.04;%delx;
C2 = 1/delx;

for i = 1:m
    for j = 1:n
           h1(i,j) = h_side(-1*mask(i,j),delx);
           h2(i,j) = h_side(mask(i,j),delx);
           b(i,j) = h1(i,j)*C1;
           if(b(i,j)~=0)
                k0(i,j) = a0;%/b(i,j);
           else
               k0(i,j) = 0;
           end
    end
end

s = signfunc(mask, absgrad, delx);

vel{1,1} = s.*h2.*gradx*C2;
vel{2,1} = s.*h2.*grady*C2;
a = h1.*k0-s.*h2.*absgrad.*cos(ca)*C2;

