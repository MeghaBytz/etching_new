function compareAnSolBiconic2D(curvconst,g,mask,pathstr,step)

% function compareAnSolBiconic2D(curvconst,g,mask,pathstr,step)

n = size(curvconst);
n = n(1)*n(2);

negative = 1;
doPlot = 0;

figure
lev = [ 0 0 ];
[ garbage, hM ] = contourf(g.xs{1}, g.xs{2}, mask, lev, 'k-'); colormap gray;
hold on;
    
C_crit = 6.67;

if(nargin >= 5 )
   %compare only one curve with analytical solution
   i = step;
   C = curvconst(i);
   if(C < C_crit)
     [data1, x] = anSolFromCurvBiconic2D(C,g,mask,negative,doPlot);
     contour(g.xs{1}, g.xs{2}, data1, lev, 'k:');
     
     fname = sprintf('%s/data_step%d',pathstr,i);
     load(fname);
     if ( ceil(i/2) == i/2 )
        contour(g.xs{1}, g.xs{2}, data, lev, 'r-');
     else
        contour(g.xs{1}, g.xs{2}, data, lev, 'g-');
     end
     axis equal
    axis(g.axis);
    hold off;
    
    figure, plotLevelSetInterior(data,0,mask);
    figure, plotLevelSetInterior(data1,0,mask);
    
    %addpath(genpath('/ices/masha/ToolboxLS-1.0/Kernel'));
    %which signedDistanceIterative
    %which addGhostExtrapolate
    %rehash
    %data = signedDistanceIterative(g,data);
    diff = abs(data-data1);
    max_abs_error = max(max(diff( find( abs(data1) < g.dx(1) ) )));
    
    diff( find( abs(data1) < g.dx(1)  )) = -1;
    diff( find( abs(data1) >=  g.dx(1) )) = 1;
    figure, plotLevelSetInterior(diff,0,mask);
    
    fprintf('\nMax abs error btw the analytical and simulated solution is %g',max_abs_error);
    
   else
     disp('No analytical solution for this step.');
   end
else
    %compare all curves that have small enough curvature
    i=2; C = curvconst(i);
    while(C < C_crit) %note - skipping step 1
        [data1, x] = anSolFromCurvBiconic2D(C,g,mask,negative,doPlot);
        contour(g.xs{1}, g.xs{2}, data1, lev, 'k:');
        if( nargin >= 4)
            fname = sprintf('%s/data_step%d',pathstr,i);
            load(fname);
            if ( ceil(i/2) == i/2 )
                contour(g.xs{1}, g.xs{2}, data, lev, 'r-');
            else
                contour(g.xs{1}, g.xs{2}, data, lev, 'g-');
            end
            
            diff = abs(data-data1);
            
            max_abs_error = max(max(diff( find( abs(data1) < g.dx(1) ) ) ));
            fprintf('\nstep %d     max_abs_error %g',i,max_abs_error);
        end
        i = i+1;
        C = curvconst(i);
    end
    
    fprintf('\n');
    axis equal
    axis(g.axis);
    hold off;    
end


