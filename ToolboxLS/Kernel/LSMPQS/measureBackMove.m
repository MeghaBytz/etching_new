function a0 = measureBackMove(pathstr,step)

% function a0 = measureBackMove(pathstr,step)

fname = sprintf('%s/curvconst.mat',pathstr);
load(fname);
fname = sprintf('%s/a.mat',pathstr);
load(fname);
fname = sprintf('%s/grid.mat',pathstr);
load(fname);
fname = sprintf('%s/mask.mat',pathstr);
load(fname);

sz = size(curvconst); sz = sz(1)*sz(2);

if (nargin == 1)
    istart = 2; istop = sz-1;
    if(istop < istart) istop=istart; end;
    
    a_back(1) = a(1);
else
    istart = step; istop = step;
end

istart
istop

outfname = sprintf('%s/measureBackMove.out',pathstr);
options = setOptions('dx',g.dx(1),'b',0.05,'OutFileName',outfname);

fid = fopen(options.OutFileName,'w');
%get the outfile directory, data will be saved in this directory
pathstr = fileparts(options.OutFileName);
if(isempty(pathstr)) 
   pathstr = '.'; %correction for current directory
end

startTime = cputime;

if( options.doDisplay )
   fprintf(fid,'Time Start %s\n',datestr(now));
   printOptions(options,fid);
end

for(i=istart:istop)
    fprintf(fid,'\nStep %d curve\n',i);
    fname = sprintf('%s/data_step%d.mat',pathstr,i);
    load(fname);
    
    backward_flag = 0;

    volume = size(find(data < 0 )); volume = volume(1);

    a0 = a(i);
    da = 0.005;
    fprintf(fid,'\nstarting curve pressure %g volume %d',a0,volume);

    while( (a0 >= 0) && (backward_flag == 0))
        a0 = a0 - da;

        options1 = setOptions(options,'a',a0,'doDisplay',0,'doSave',0, 'doReinit',1);
        fclose(fid);
        [data1, g, mask, curv1] = constCurvModel(options1,data,g,mask);
        fid = fopen(options.OutFileName,'a'); 

        volume1 = size(find(data1 < 0 )); volume1 = volume1(1);

        if( volume1 < volume )
                backward_flag = 1;
                fprintf(fid,'\na %g volume1 %d backward_flag %d\n',a0,volume1,backward_flag);
                if( istop > istart )
                    a_back(i) = a0;
                end;
        end

        %fprintf(fid,'\na %g volume1 %d backward_flag %d\n',a0,volume1,backward_flag);
    end
end

if( istop > istart )
    fname = sprintf('%s/a_back.mat',pathstr);
    save(fname,'a_back');

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
            H2 = plot(k_avg,diff,'r.'), xlabel('curvature'); ylabel('rel. press. diff.');
            legend([ H1(1), H2(1)], 'constcurv', 'computed kavg');
end

fclose(fid);