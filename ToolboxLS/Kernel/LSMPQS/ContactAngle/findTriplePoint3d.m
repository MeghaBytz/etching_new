function point = findTriplePoint(data,g,mask)
% data - level set function, interested in level = 0
% g - grid
% mask - level set function defining pore (and grain) space

 sz1 = size( find(data < 0)); sz1 = sz1(1);
 sz2 = size( find((data > 0) & (mask > 0)) ); sz2 = sz2(1);
 sz3 = size( find((data > 0) & (mask < 0)) ); sz3 = sz3(1);

 point = 0;
 
 if ( sz1*sz2*sz3 == 0 )
     disp('findTriplePoint() - data and mask do not define 3 phases. Abort');
     return;
 end
 
 C = contour(data,[0 0]); % find a 0-level set and output coords to matrix C
 Cmask = contour(mask,[0 0]);
 
 C1 = getContourPoints(C); C1 = C1';
 Cmask1 = getContourPoints(Cmask); Cmask1 = Cmask1';
 
 [grain_intfc i1 i2] = intersect(C1,Cmask1,'rows');
 x_grn = grain_intfc(:,1);
 y_grn = grain_intfc(:,2);
 
 % fluid interface are points on boundary of data, that are not on
 % grain(mask) boundary
 fld_intfc = setdiff(C1,grain_intfc,'rows');
 x_fld = fld_intfc(:,1);
 y_fld = fld_intfc(:,2);
 
 %[x_edge1 y_edge1 count1] = findEdgePoint(grain_intfc,data,mask);
 [x_edge2 y_edge2 count2] = findEdgePoint(fld_intfc,data,mask);
 
 [xc1 yc1] = findClosestPoint(grain_intfc,data,mask,x_edge2,y_edge2);
 [xm1 ym1] = findClosestPoint(fld_intfc,data,mask,xc1,yc1);
 
 yc1 = yc1';
 ym1 = ym1';
 
 grad_data_interp = interp2(grad_norm_data,xc1,yc1)
 grad_mask_interp = interp2(grad_norm_mask,xc1,yc1)
 dx_data_interp   = interp2(    dx_data,xc1,yc1);
 dy_data_interp   = interp2(    dy_data,xc1,yc1);
 dx_mask_interp   = interp2(    dx_mask,xc1,yc1);
 dy_mask_interp   = interp2(    dy_mask,xc1,yc1);
 
 
 cos_angle = (dx_data_interp .* dx_mask_interp + dy_data_interp .* dy_mask_interp );
 cos_angle  = cos_angle ./(grad_data_interp .* grad_mask_interp)
 

 if(g.dim == 2)
     %plot points
     figure
     plot(y_grn,x_grn,'.'),xlabel('x'), ylabel('y'), hold on
     plot(y_fld,x_fld,'r.'), hold on
     
     plot(yc1,xc1,'m.'), hold on
     plot(ym1,xm1,'k.'), hold on
    
     if( count2 > 0)
         plot(y_edge2,x_edge2,'g.')
     end
     hold off
     axis([0 g.N(1) 0 g.N(2)])
 end
 
 
 
 
 function [x_edge y_edge count] = findEdgePoint(C,data,mask)
 %find point that has 3 phases in the neighborhood
 %works in 2D for now
 
 [N tmp] = size(C);
 [m n] = size(data);
 
 count = 0; x_edge(1) = 0; y_edge(1) = 0;
 for(i=1:N)
    x = C(i,1); y = C(i,2);
    k = floor(x); l = floor(y);
    %fprintf('\ni= %d x %g y %g k %d l %d',i,x,y,k,l);
    [grain inside outside] = test_phase(data,mask,l,k);
    
    if( k+1 <= n)
        [gr in out] = test_phase(data,mask,l,k+1);
        grain = grain | gr; inside = inside | in; outside = outside | out;
    end;
    
    if( l+1 <= m )
        [gr in out] = test_phase(data,mask,l+1,k);
        grain = grain | gr; inside = inside | in; outside = outside | out;
        
        if( k+1 <= n)
             [gr in out] = test_phase(data,mask,l+1,k+1);
             grain = grain | gr; inside = inside | in; outside = outside | out;
        end;
    end;
    
    if( grain & inside & outside )
        count = count+1;
        x_edge(count) = x;
        y_edge(count) = y;
    end
 end
 
 
 function [grain inside outside] = test_phase(data,mask,i,j)
 % test which image phase does the point belong to
 grain = 0; inside = 0; outside = 0;
 small = 10;
 
 if( mask(i,j) > -small*eps ) grain = 1;
 elseif( data(i,j) > -small*eps) outside = 1;
     %if( mask(i,j) > 0 ) grain = 1;
 %elseif( data(i,j) > 0 ) outside = 1;
 else inside = 1;
 end
 
 %fprintf('\n\tmask %g data %g grain %d inside %d outside %d',mask(i,j),data(i,j),grain,inside,outside);
 
 function [xc yc] = findClosestPoint(C,data,mask,x,y)
 %find point n C closest to x,y
 %works in 2D for now
 
 [N tmp] = size(C);
 [m n] = size(x); len = m*n; % either m or n is 1
 
 for(j=1:len)
    xc(j) = 0; yc(j) = 0; dist = 1000;
    for(i=1:N)
        xp = C(i,1); yp = C(i,2);

        dist1 = (x(j) - xp)*(x(j)-xp) + (y(j)-yp)*(y(j)-yp);
        if dist1 < dist
            xc(j) = xp;
            yc(j) = yp;
            dist = dist1;
        end
    end
 end
 