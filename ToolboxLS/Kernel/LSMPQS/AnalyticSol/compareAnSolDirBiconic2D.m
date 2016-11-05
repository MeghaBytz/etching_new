function compareAnSolDirBiconic2D(pathstr,step)

% compareAnSolDirBiconic2D(pathstr,step)
% compareAnSolDirBiconic2D(pathstr)

fname = sprintf('%s/curvconst.mat',pathstr);
load(fname);
fname = sprintf('%s/grid.mat',pathstr);
load(fname);
fname = sprintf('%s/mask.mat',pathstr);
load(fname);

if( nargin >= 2 )
    compareAnSolBiconic2D(curvconst,g,mask,pathstr,step);
else
    compareAnSolBiconic2D(curvconst,g,mask,pathstr);
end