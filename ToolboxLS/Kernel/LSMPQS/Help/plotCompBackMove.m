function [rel_press abs_press rel_press_deriv] = plotCompBackMove(pathstr,crit_step,num)

% function replotBackMove(pathstr)

fname = sprintf('%s/curvconst.mat',pathstr);
load(fname);
fname = sprintf('%s/a.mat',pathstr);
load(fname);

sz = size(curvconst); sz = sz(1)*sz(2);

istart = 2; istop = sz-1;
if(istop < istart) istop=istart; end;
a_back(1) = a(1);

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
    
    figure, H1=plot(i,rel_diff,'b.'), xlabel('Step'); ylabel('rel. press. diff.');
            hold on;
            
    for(j=1:num)  
      j1 = crit_step(j) - istart + 1;
      plot(crit_step(j),rel_diff(j1),'m*');
      rel_press(j) = rel_diff(j1);
      abs_press(j) = rel_diff(j1)*a(j1);
      
      step_diff = 2;
      if( j1+1 <= numel(rel_diff) ) 
          rel_press_deriv(j) = rel_diff(j1+1);       
      else
          rel_press_deriv(j) = rel_diff(j1);
          step_diff = step_diff - 1;
      end
      
      if( j1-1 >= 1 )
          rel_press_deriv(j) = rel_press_deriv(j) - rel_diff(j1-1);       
      else
          rel_press_deriv(j) = rel_press_deriv(j) - rel_diff(j1);
          step_diff = step_diff - 1;
      end
      
      if(step_diff) 
          rel_press_deriv(j) = rel_press_deriv(j) / step_diff;
      end
      
    end
    
else error('replotBackMove() - Nothing to plot.');
end
