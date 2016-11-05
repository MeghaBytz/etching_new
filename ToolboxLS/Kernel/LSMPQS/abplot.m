function [a,b,vel] = abplot(mask,delx,a0,ca)
h = @heaviside;
s = @signfunc;

[ny,nx] = size(mask);
a(1:ny,1:nx) = 0;
b(1:ny,1:nx) = h(-1*mask(1:ny,1:nx));
vel{2,1}(ny,nx) = 0;
vel{1,1}(ny,nx) = 0;
k0 = a0./b;
%
mask1(2:ny+1,2:nx+1)=mask;
mask1(1,2:nx+1)=mask(1,:);
mask1(ny+2,2:nx+1)=mask(ny,:);
mask1(2:ny+1,1)=mask(:,1);
mask1(2:ny+1,nx+2)=mask(:,nx);

C = 1/delx;

%
for i=2:ny+1
    for j=2:nx+1
        [absgr, gradx, grady] = gradmask(i,j);          %recalculating a,b
        a(i-1,j-1) = h(-1*mask1(i,j))*k0(i-1,j-1)-C*s(mask1(i,j),absgr)*h(mask1(i,j))*cos(pi-ca)*absgr;
        vel{1,1}(i-1,j-1) = C*s(mask1(i,j),absgr)*h(mask1(i,j))*gradx;
        vel{2,1}(i-1,j-1) = C*s(mask1(i,j),absgr)*h(mask1(i,j))*grady;
    end
end


    function y = heaviside(x)
        eps = 1.5*delx;
        if (x<(-1*eps))
            y = 0;
        elseif (x>(eps))
            y = 1;
        else
            y = 0.5 + x/(2*eps) + (0.5/pi)*sin((pi*x)/(eps));
        end
    end

    function [absgrad,gradx,grady]= gradmask(u,v)
        gradx = 0.5*(mask1(u+1,v)-mask1(u-1,v))/delx;           %changed (u,v+1) to (u+1,v) etc., as grid is not as direction.
        grady = 0.5*(mask1(u,v+1)-mask1(u,v-1))/delx;
        absgrad = (gradx^2+grady^2)^0.5;
    end

    function k = signfunc(w,absgr)
        k = w/(w^2+(absgr*delx)^2)^0.5;
    end
end



%}