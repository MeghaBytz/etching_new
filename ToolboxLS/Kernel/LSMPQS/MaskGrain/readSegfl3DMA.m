function [data nx ny nz] = readSegfl3DMA(filename)

% [data nx ny nz] = readSegfl3DMA(filename)
% filename - basename or full name of 3DMA segmented file to be read
% reads 3DMA segmented (volume) file into a 3D array 'data' whose size is
% [nx ny nz] and values are 0 and 1, this array is double by default

do_vis = 0; % set to one if quick visualization of isosurface desired

if( findstr('.gz',filename) )
    zip = 1;
    filename = gunzip(filename);
    filename1 = filename{1};
else
    filename1 = filename;
    zip = 0;
end

fid = fopen(filename1,'r');
e  = fread(fid,1,'char');
nx = fread(fid,1,'int');
ny = fread(fid,1,'int');
zs = fread(fid,1,'int');
ze = fread(fid,1,'int');
nz = (ze - zs + 1);
bnxyz = fread(fid,1,'int');
nresid = fread(fid,1,'int'); 

bitpacked = fread(fid,bnxyz+1,'uint8');
data = [];
% unpack into a row array 'data'
for(i=1:(bnxyz+1))
    a = bitget(bitpacked(i),1:8); 
    %note that bitget returns bits in reverse order, need to flip it
    a = fliplr(a);
    data = [data, a];   
end

clear bitpacked;


nxy = nx*ny; nxyz = nxy*nz;
%only the first nxyz elements from bitpacked array are relevant
data = data(1:nxyz);
data = reshape(data,[nx ny nz]);

%temporary, want things turned around
data = data';
[nx ny nz] = size(data);

%testing and quick visualization
if( do_vis & nz > 1)
    % 3D visualization
    x = [1:nx]; y = [1:ny]; z = [1:nz];
    [X Y Z] = meshgrid(x,y,z);
    figure, h = patch(isosurface(X,Y,Z,data,0.5));
    set(h,'FaceColor','blue','EdgeColor', 'none');
    daspect([1 1 1]); camlight; lighting phong
    view(3)
elseif (do_vis)
    % 2D visualization
    figure, contourf(data,[0.5 0.5],'k-'); colormap gray;
end

fclose(fid);

if( zip )
    cmnd = sprintf('/bin/rm %s',filename1);
    system(cmnd);
end
