function analyzeCurvature1(pathstr,fid)

%function analyzeCurvature1(pathstr)
%function analyzeCurvature1(pathstr,fid)

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

% count voxels in connected components of the invading fluid at entry
dim = 1;
[m n] = size(mask);
min_cross_sec = m;
for(i=1:m)
    [ind comp sz_comp] = findConnComp(mask,dim,i);
    cross_sec = min(sz_comp)*g.dx(1);
    if(cross_sec < min_cross_sec)
        min_cross_sec = cross_sec;
    end
end
radius = min_cross_sec/2.0;
Cexact = 1/(radius);

%real geometries might need a fix.
%Cexact = 2/sqrt(0.06*0.06 + 0.12*0.12);
%Cexact = 2/sqrt(0.6*0.6 + 0.2*0.2);
%Cexact = 1/0.1;

%Cexact = 1/0.08485;
Cexact = 1/0.075
%The first (initializing) step is skipped
n = size(k_max);
n = n(1)*n(2);
step = (2:(n-1));

crit_step = n-1;
found = 1; i = crit_step;
while(found)
    if( k_avg(i) <= k_avg(i-1))
        i= i-1;
    else
        found = 0;
    end
end
crit_step = i;

plotCompCurv1(k_max(2:(n-1)),k_avg(2:(n-1)),curvconst(2:(n-1)),step,Cexact,crit_step-1);
plotCompData(pathstr,1,n,crit_step);

crit_curv = k_avg(crit_step);
abs_error = abs(Cexact - crit_curv);
rel_error = abs_error/Cexact;

crit_curv1 = curvconst(crit_step);
abs_error1 = abs(Cexact - crit_curv1);
rel_error1 = abs_error1/Cexact;


fprintf('\n%s',pathstr);
fprintf('\nexact crit curv. %g',Cexact);
fprintf('\ncrit_step k_avg abs_error rel_error');
fprintf('\n%d        %2.6g  %2.6g  %2.6g', crit_step, crit_curv, abs_error, rel_error);
fprintf('\n          k_simul abs_error rel_error');
fprintf('\n          %2.6g  %2.6g  %2.6g', crit_curv1, abs_error1, rel_error1);
    
if(nargin >= 2)
    fprintf(fid,'\n%s',pathstr);
    fprintf(fid,'\nexact crit curv. %g',Cexact);
    fprintf(fid,'\ncrit_step k_avg abs_error rel_error');
    fprintf(fid,'\n%d        %2.6g  %2.6g  %2.6g', crit_step, crit_curv, abs_error, rel_error);
    fprintf(fid,'\n          k_simul abs_error rel_error');
    fprintf(fid,'\n          %2.6g  %2.6g  %2.6g', crit_curv1, abs_error1, rel_error1);
else
    fname = sprintf('%s/error',pathstr);
    fid = fopen(fname,'w');
    fprintf(fid,'\n%s',pathstr);
    fprintf(fid,'\nexact crit curv. %g',Cexact);
    fprintf(fid,'\ncrit_step k_avg abs_error rel_error');
    fprintf(fid,'\n%d        %2.6g  %2.6g  %2.6g', crit_step, crit_curv, abs_error, rel_error);
    fprintf(fid,'\n          k_simul abs_error rel_error');
    fprintf(fid,'\n          %2.6g  %2.6g  %2.6g', crit_curv1, abs_error1, rel_error1);
    fclose(fid);
end

fprintf('\n');
