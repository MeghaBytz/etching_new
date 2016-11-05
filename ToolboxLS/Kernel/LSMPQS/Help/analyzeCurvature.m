function analyzeCurvature(pathstr,fid)

%function analyzeCurvature(pathstr)
%function analyzeCurvature(pathstr,fid)

fname = sprintf('%s/curvconst.mat',pathstr);
load(fname);
fname = sprintf('%s/k_avg.mat',pathstr);
load(fname);
fname = sprintf('%s/k_max.mat',pathstr);
load(fname);

fname = sprintf('%s/mask.mat',pathstr);
load(fname);
fname = sprintf('%s/grid.mat',pathstr);
load(fname);


%for(i=1:g.N(2))
%    tmp = size(find(mask(i,:) < 0));
%    cross(i) = tmp(1) *tmp(2);
%end
%min_cross_sec = min(cross)*g.dx(1);
%radius = min_cross_sec/2.0;

%Cexact = 1/(radius);
Cexact = 1/0.15;
%Cexact = 2/sqrt(0.06*0.06 + 0.12*0.12);
%Cexact = 1/0.07;
%Cexact = 1/0.1;

%The first (initialzing) step is skipped
n = size(k_max);
n = n(1)*n(2)
step = (2:n);
%Cexact=2*Cexact;
[crit_curv crit_step abs_error rel_error num] = plotCompCurv(k_max(2:n),k_avg(2:n),curvconst(2:n),step,Cexact);
plotCompData(pathstr,n,crit_step);

fname = sprintf('%s/a_back.mat',pathstr);
if( exist(fname,'file'))
    [rel_press abs_press rel_press_deriv] = plotCompBackMove(pathstr,crit_step,num);
else
    rel_press = zeros(size(crit_step));
    abs_press = zeros(size(crit_step));
    rel_press_deriv = zeros(size(crit_step));
end

fprintf('\n%s',pathstr);
fprintf('\nexact crit curv. %g',Cexact);
fprintf('\ncrit_step crit_curv abs_error rel_error abs_press rel_press rel_press_deriv');
for(i=1:num)
    fprintf('\n%8d  %2.6g  %2.6g  %2.6g  %2.6g  %2.6g  %2.6g', crit_step(i), crit_curv(i),...
                        abs_error(i), rel_error(i), rel_press(i), abs_press(i), rel_press_deriv(i));
end
    
if(nargin > 1)
    fprintf(fid,'\n%s',pathstr);
    %fprintf(fid,'\nexact crit curv. %g',Cexact);
    fprintf(fid,'\ncrit_step   crit_curv  abs_error rel_error  abs_press   rel_press   rel_press_deriv');
    for(i=1:num)
        fprintf(fid,'\n%8d  %2.6g  %2.6g  %2.6g  %2.6g  %2.6g  %2.6g',crit_step(i), crit_curv(i),...
                            abs_error(i), rel_error(i), rel_press(i), abs_press(i), rel_press_deriv(i));
    end
end

fprintf('\n');
