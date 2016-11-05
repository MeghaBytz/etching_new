function [crit_curv crit_step abs_error rel_error num] = plotCompCurv(k_max,k_avg,C,step,Cexact)

% function plotCompCurv(k_max,k_avg,C,step,Cexact)
% k_avg - average curvature array
% k_max - max curvature array
% C - simulation curvature array
% step - step id's to put on x-axis
% Cexact (scalar) - (optional) exact critical curvature

figure, h1 = plot(step,k_max, 'b.'); hold on
xlabel('Step'), ylabel('Curvature at level set tip')
h2 = plot(step,k_avg, 'r.');
h3 = plot(step,C,'k.')

if(nargin >= 4)
    const = Cexact * ones(size(k_max));
    h4 = plot(step,const,'g');
    legend([ h1(1), h2(1), h3(1),h4(1) ], ...
         {'kmax', 'kavg','LSM curvature','critical'});
else
     legend([ h1(1), h2(1), h3(1)], ...
         {'kmax', 'kavg','LSM curvature'});
end

n = size(k_avg);
n = n(1)*n(2);

step_s = 0; step_e = 0; found = 0;

i=2;
num = 0;
while(i<=n-1)
    if( ( ( k_avg(i-1) > k_avg(i) ) && (k_avg(i) >= k_avg(i+1) ) ||...
        ( k_avg(i-1) >= k_avg(i) ) && (k_avg(i) > k_avg(i+1) ))  && (k_avg(i+1) > 0) )
       step_s = i-1;
       step_e = i+1;
       j=i+2;
       while( (j<=n) && (k_avg(j) < k_avg(j-1))  )
          step_e = j;
          j = j+1;
       end
       
       num = num + 1;
       crit_curv(num) = k_avg(step_s);
       crit_step(num) = step(step_s);
       abs_error(num) = abs(k_avg(step_s)-Cexact);
       rel_error(num) = abs_error(num)/Cexact;
       fprintf('critical curvature value k_avg(%d) = %g rel.error %g\n',step(step_s),k_avg(step_s),rel_error(num));
       plot(step(step_s:step_e),k_avg(step_s:step_e),'m--');
       i = step_e;
    end
    i = i+1;
end

if( num == 0 )
    % well, search again with lower criteria
    disp('lower criteria search');
    i=2;
    while(i<=n-1)
        if( ( k_avg(i-1) > k_avg(i) ) && (k_avg(i) > 0) )
           step_s = i-1;
           step_e = i+1;
           j=i+2;
           while( (j<=n) && (k_avg(j) < k_avg(j-1))  )
              step_e = j;
              j = j+1;
           end

           num = num + 1;
           crit_curv(num) = k_avg(step_s);
           crit_step(num) = step(step_s);
           abs_error(num) = abs(k_avg(step_s)-Cexact);
           rel_error(num) = abs_error(num)/Cexact;
           fprintf('critical curvature value k_avg(%d) = %g rel.error %g\n',step(step_s),k_avg(step_s),rel_error(num));
           plot(step(step_s:step_e),k_avg(step_s:step_e),'m--');
           i = step_e;
        end
    i = i+1;
end

  
end    