function plotCompDataGray(pathstr,start,num,crit_step,skip_step)

%function plotCompDataGray(pathstr,start,num,crit_step,skip_step)
% replots data, start = starting step, num = ending step
% standard filenames assumed (see below)
% crit_step -- steps to be plot in thicker line
% skip_step - step to be skipped

fname = sprintf('%s/mask.mat',pathstr);
load(fname);
fname = sprintf('%s/grid.mat',pathstr);
load(fname);
fname = sprintf('%s/data_init.mat',pathstr);
load(fname);

%set flip to 1 if you want flip x and y axis on the plot
flip = 1;
if( flip )
    mask = mask';
    data0 = data0';
    tmp = g.xs{1}; g.xs{1} = g.xs{2}; g.xs{2} = tmp;
    g.xs{1} = g.xs{1}';
    g.xs{2} = g.xs{2}';
    tmp1 = g.axis;
    g.axis(1) = g.axis(3);
    g.axis(2) = g.axis(4);
    g.axis(3) = tmp1(1);
    g.axis(4) = tmp1(2);
    tmp2 = g.N(1); g.N(1 ) = g.N(2); g.N(2) = tmp2;
    tmp2 = g.shape(1); g.shape(1) = g.shape(2);  g.shape(2) = tmp2; 
end

if (nargin < 4)
    crit_step = 0;
end

if (nargin < 5)
    skip_step = 0;
end

figure;
  lev = [ 0 0 ];
  if(exist('data'))
    [ garbage, hI ] = contour(g.xs{1}, g.xs{2}, data, lev,'k:');
  elseif(exist('data0'))
    [ garbage, hI ] = contour(g.xs{1}, g.xs{2}, data0, lev,'k:');
  end
  hold on;
  [ garbage, hM ] = contourf(g.xs{1}, g.xs{2}, mask, lev);colormap gray;
  
  hs = [ hI];
  %set(hs, 'LineWidth', 2.0);
  
  axis equal
  axis(g.axis);

j = 1;
for(i=start:num)
    if( i ~= skip_step )
        fname = sprintf('%s/data_step%d.mat',pathstr,i);
        load(fname);
        
        if (flip ) data = data'; end;

        if((j <= numel(crit_step)) && (i == crit_step(j)))
            [ garbage, hF ] = contour(g.xs{1}, g.xs{2}, data, lev, 'k-', 'LineWidth', 2.0); hold on;
            j = j+1;
            legend([ hI(1), hF(1) ], {'initial','critical'}, 2);
        elseif ( ceil(i/2) == i/2 )
            [ garbage, hF ] = contour(g.xs{1}, g.xs{2}, data, lev, '-','Color', [0.6 0.6 0.6]); hold on;
        else
            [ garbage, hF ] = contour(g.xs{1}, g.xs{2}, data, lev, 'k-'); hold on;
        end

        %set(hF, 'LineWidth', 2.0);
        axis equal
        axis(g.axis);
    end
end
hold off
    
    
    
    
    


