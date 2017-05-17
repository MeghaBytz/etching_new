function plotKavg(k_max,k_avg,C,Cexact)

% function plotKavg(k_max,k_avg,C)
% function plotKavg(k_max,k_avg,C,Cexact)
% k_avg - average curvature array
% k_max - max curvature array
% C - simulation curvature array
% Cexact (scalar) - (optional) exact critical curvature

figure, h1 = plot(k_max, 'b.'); hold on
xlabel('Step'), ylabel('Curvature at level set tip')
h2 = plot(k_avg, 'r.');
h3 = plot(C,'k.');

if(nargin >= 4)
    const = Cexact * ones(size(k_max));
    h4 = plot(const,'g');
    legend([ h1(1), h2(1), h3(1),h4(1) ], ...
         {'kmax', 'kavg','LSM curvature','critical'});
else
     legend([ h1(1), h2(1), h3(1)], ...
         {'kmax', 'kavg','LSM curvature'});
end