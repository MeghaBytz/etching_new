function compareAnSolDirThroat2D(pathstr,step)

% function compareAnSolDirThroat2D(pathstr,step)
% function compareAnSolDirThroat2D(pathstr)

fname = sprintf('%s/curvconst.mat',pathstr);
load(fname);
fname = sprintf('%s/grid.mat',pathstr);
load(fname);
fname = sprintf('%s/mask.mat',pathstr);
load(fname);

if( nargin >= 2 )
    compareAnSolThroat2D(curvconst,g,mask,pathstr,step);
else
    compareAnSolThroat2D(curvconst,g,mask,pathstr);
end