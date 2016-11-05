function replotBackMove(pathstr)

% function replotBackMove(pathstr)

fname = sprintf('%s/curvconst.mat',pathstr);
load(fname);
fname = sprintf('%s/a.mat',pathstr);
load(fname);

sz = size(curvconst); sz = sz(1)*sz(2);

if (nargin == 1)
    istart = 2; istop = sz-1;
    if(istop < istart) istop=istart; end;
    
    a_back(1) = a(1);
else
    istart = step; istop = step;
end
i = (istart:istop);

if( istop > istart )
    fname = sprintf('%s/a_back.mat',pathstr);
    load(fname);

    a = a(istart:istop);
    a_back = a_back(istart:istop);
    diff = a-a_back;
    rel_diff = abs( (a-a_back)./a );
    curvconst = curvconst(istart:istop);

    fname = sprintf('%s/k_avg.mat',pathstr);
    load(fname);
    k_avg = k_avg(istart:istop);
    
    figure, H1= plot(curvconst,diff,'b.'), xlabel('curvature'); ylabel('pressure diff.');
            hold on;
            H2 = plot(k_avg,diff,'r.'), xlabel('curvature'); ylabel('pressure diff.');
            legend([ H1(1), H2(1)], 'constcurv', 'computed kavg');
    figure, H1=plot(curvconst,rel_diff,'b.'), xlabel('curvature'); ylabel('rel. press. diff.');
            hold on;
            H2 = plot(k_avg,rel_diff,'r.'), xlabel('curvature'); ylabel('rel. press. diff.');
            legend([ H1(1), H2(1)], 'constcurv', 'computed kavg');
    figure, H1=plot(i,rel_diff,'b.'), xlabel('Step'); ylabel('rel. press. diff.');
            hold off;
    
else error('replotBackMove() - Nothing to plot.');
end
