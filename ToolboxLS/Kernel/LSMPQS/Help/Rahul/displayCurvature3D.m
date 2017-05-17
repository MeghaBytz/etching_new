function [k_avg k_min k_max] = displayCurvature3D(data,g,mask,curv,fid,plot)

% function [k_avg k_min k_max] = displayCurvature(data,g,mask,curv,fid,plot)
% data - level set function, will plot curvature values in 'curv' on its level = 0
% g    - grid
% mask - level set function defining pore and grain space
% curv - curvature data
% fid  - (optional) file id of an outfile
%       Outfile, if supplied, will be appended with curvature average, min, max.
% plot - (optional) 1 if you want a figure with curvature plot, defaults to 1

%compute current curvature
 sz = size( find(data < 0)); sz = sz(1);
 sz1 = size( find(data > 0)); sz1 = sz1(1);
 
 if (nargin < 6)
     plot = 1; %plot by default
 end

 if ( sz *sz1 == 0 )
     k_avg = 0; k_min = 0; k_max = 0;
     return;
 end
 
 C = contourc(data,[0 0]); % find a 0-level set and output coords to matrix C
 Cmask = contourc(mask,[0 0]);
 
 C1 = getContourPoints(C); C1 = C1';
 Cmask1 = getContourPoints(Cmask); Cmask1 = Cmask1';  
 
 [common i1 i2] = intersect(C1,Cmask1,'rows');
 % fluid interface are points on boundary of data, that are not on
 % grain(mask) boundary
 fld_intfc = setdiff(C1,common,'rows');
 
 x_fld = fld_intfc(:,1);
 y_fld = fld_intfc(:,2);
 y_fld = y_fld';
 
 %find interpolated curvature values at these coords
 curv_interp = interp2(curv,x_fld,y_fld); % this creates interpolated values on entire grid defined by x,y
 K = diag(curv_interp); %isolate only values at (x_fld(i),y_fld(i))
 sz = size(K); sz = sz(1);
 
 if( sz )
     k_avg_fld_intfc = mean(K);
     k_min_fld_intfc = min(K);
     k_max_fld_intfc = max(K);
 else
   k_avg_fld_intfc  = 0; k_min_fld_intfc  = 0; k_max_fld_intfc  = 0;
 end
   
 % this works for straight boundaries - tube
 %[M N] = size(data);
 %y_fld = y_fld';
 %I = find( (x_fld >= 5) & (x_fld <= N-5) & (y_fld >= 5) & (y_fld <= M-5) ); %5 grid points away from volume sides; works well for duct
 %y_fld = y_fld';
 %K1 = K(I);
 
 % points close to pore space boundary should be avoided
 dist = -2*g.dx(1);
 mask_interp = interp2(mask,x_fld,y_fld);
 mask_interp = diag(mask_interp);
 %min1 = min( mask_interp )
 %max1 = max( mask_interp )
 I = find( mask_interp < dist );
 K1 = K(I); x_fld1 = x_fld(I); y_fld1 = y_fld(I);
 
 sz1 = size(K1); sz1 = sz1(1);
 
 if( sz1 ) 
     k_avg = mean(K1) ;
     k_min = min(K1) ;
     k_max = max(K1);
 else
     k_avg = 0; k_min = 0; k_max = 0;
 end
 
 if (nargin >= 5)
     fprintf(fid,'Fluid interface curvature ( data = 0) & (mask < 0) - %d values\n',sz);
     fprintf(fid,'\tk_avg_fld_intfc %g\n',k_avg_fld_intfc);
     fprintf(fid,'\tk_min_fld_intfc %g\n',k_min_fld_intfc);
     fprintf(fid,'\tk_max_fld_intfc %g\n',k_max_fld_intfc);
      fprintf(fid,'Fluid interface curvature values 2 grid pts away from volume bdry - %d values\n',sz1);
     fprintf(fid,'\tk_avg %g\n',k_avg);
     fprintf(fid,'\tk_min %g\n',k_min);
     fprintf(fid,'\tk_max %g\n',k_max);
 else
     fprintf('Fluid interface curvature ( data = 0) & (mask < 0) - %d values\n',sz);
     fprintf('\tk_avg_fld_intfc %g\n',k_avg_fld_intfc);
     fprintf('\tk_min_fld_intfc %g\n',k_min_fld_intfc);
     fprintf('\tk_max_fld_intfc %g\n',k_max_fld_intfc);
      fprintf('Fluid interface curvature values 3 grid pts away from volume bdry - %d values\n',sz1);
     fprintf('\tk_avg %g\n',k_avg);
     fprintf('\tk_min %g\n',k_min);
     fprintf('\tk_max %g\n',k_max);
 end
 
 if(g.dim == 2 && plot > 0)
     %plot curvature  
     figure,plot3(y_fld,x_fld,K,'.'),xlabel('x'), ylabel('y'), zlabel('fld intfc curvature')
     hold on
     plot3(y_fld1,x_fld1,K1,'r.'),
     axis([0 g.N(1) 0 g.N(2) k_min_fld_intfc k_max_fld_intfc])
 end
