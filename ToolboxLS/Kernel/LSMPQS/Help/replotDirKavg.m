function replotDirKavg(pathstr)

% function replotDirKavg(pathstr)

fname = sprintf('%s/curvconst.mat',pathstr);
load(fname);
fname = sprintf('%s/k_avg.mat',pathstr);
load(fname);
fname = sprintf('%s/k_max.mat',pathstr);
load(fname);

plotKavg(k_max,k_avg,curvconst);