function readMergePlots(pathstr,step1,step2,step3,step4,step5)

% function readMergePlots(pathstr,step1,step2,step3,step4,step5)

fname = sprintf('%s/mask.mat',pathstr);
load(fname);
fname = sprintf('%s/grid.mat',pathstr);
load(fname);

fname = sprintf('%s/data_init.mat',pathstr);
load(fname);
fname = sprintf('%s/curvconst.mat',pathstr);
load(fname);

% temporary scaling
scale_fac = 1.0;
g.dx = g.dx/scale_fac;
g.min = g.min/scale_fac;
g.max = g.max/scale_fac;
g.xs{1} = g.xs{1}/scale_fac;
g.xs{2} = g.xs{2}/scale_fac;
g.vs{1} = g.vs{1}/scale_fac;
g.vs{2} = g.vs{2}/scale_fac;
g.axis = g.axis/scale_fac;
curvconst = curvconst*scale_fac;

fname = sprintf('%s/data_step%d.mat',pathstr,step1);
load(fname); data1 = max(mask,data);
%d1_title = sprintf('step %d',step1);
d1_title = sprintf('c = %2.2f',curvconst(step1));

fname = sprintf('%s/data_step%d.mat',pathstr,step2);
load(fname); data2 = max(mask,data);
%d2_title = sprintf('step %d',step2);
d2_title = sprintf('c = %2.2f',curvconst(step2));

fname = sprintf('%s/continue/data_step%d.mat',pathstr,step3);
load(fname); data3 = max(mask,data);
%d3_title = sprintf('step %d',step3);
d3_title = sprintf('c = %2.2f',0.0391);
%d3_title = sprintf('c = %2.2f',curvconst(step3));
if(nargin <= 4)
    % d3 will have name "final"
    mergePlots(data0,data1,data2,data3,g,mask,d1_title,d2_title,d3_title);
else
    fname = sprintf('%s/data_step%d.mat',pathstr,step4);
    load(fname); data4 = max(mask,data);
    %d4_title =  sprintf('step %d',step4);
    d4_title = sprintf('c = %2.2f',curvconst(step4));
    
    fname = sprintf('%s/data_step%d.mat',pathstr,step5);
    load(fname); data5 =max(mask,data);
    %d5_title =  sprintf('step %d',step5);
    d5_title = sprintf('c = %2.2f',curvconst(step5));
    
    %data4 will be called "critical'" and data5 will be called "final"
    mergePlots1(data0,data1,data2,data3,data4,data5,g,mask,d1_title,d2_title,d3_title,d4_title,d5_title);
end
    









