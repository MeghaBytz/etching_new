function plotCompCurv1(k_max,k_avg,C,step,Cexact,crit_step)

% function plotCompCurv1(k_max,k_avg,C,step,Cexact,crit_step)
% k_avg - average curvature array
% k_max - max curvature array
% C - simulation curvature array
% step - step id's to put on x-axis
% Cexact (scalar) - (optional) exact critical curvature

figure, h1 = plot(step,k_max, 'b.'); hold on
xlabel('Step'), ylabel('Curvature at level set tip')
h2 = plot(step,k_avg, 'r.');
h3 = plot(step,C,'k.')
if(nargin >= 6)
    h4 = plot(step(crit_step),k_avg(crit_step),'m*');
end

if(nargin >= 4)
    const = Cexact * ones(size(k_max));
    h4 = plot(step,const,'g');
    legend([ h1(1), h2(1), h3(1),h4(1) ], ...
         {'kmax', 'kavg','LSM curvature','critical'});
else
     legend([ h1(1), h2(1), h3(1)], ...
         {'kmax', 'kavg','LSM curvature'});
end

