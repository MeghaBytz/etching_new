function [data, data0, g, mask] = initializeDataGrid(options,MaskGeom,periodic)

% function [data, data0, g, mask] = initializeDataGrid(options,MaskGeom,periodic)
% initializes data,data0, grid, mask
% options - structure with input options, see setOptions.m
% MaskGeom - list returned from setAvailGeom.m
% periodic - LSToolbox supports periodic boundary condition, but hasn't been tested; assumed 0

dxIn = options.dx;
testGeom = options.testGeom;
ICplane = options.ICplane;
testGeom1 = options.testGeom1;
testGeom2 = options.testGeom2;   testGeom3 = options.testGeom3;
patchx = options.PatchX;
patchy = options.PatchY;

nx = ceil(1/dxIn); % need to avoid rounding off errors
g.dx = 1/nx %input value

% specify basic grid info for each geometry
if( testGeom == 1)
    % see maskDuct2D.m
    g.dim = 2;
    g.min = [-1 - g.dx; 0 - g.dx];
    g.max = +1 + g.dx;
elseif( testGeom == 2)
    %(see maskCircDuct2D.m
    g.dim = 2;
    g.min = [-0.5 - g.dx; -0.2 - g.dx];
    g.max = [ 0.5 + g.dx;  0.3 + g.dx];
elseif( testGeom == 3)
    %maskSphere2D.m
    g.dim = 2;
    g.min = [0-g.dx;0 - g.dx];
    g.max = [1 + g.dx;1.16 + g.dx]; 
elseif( testGeom == 28   )
    %
    g.dim = 2;
    g.min = [-1.0;-0.5];
    g.max = [2.76 + g.dx;2.7 + g.dx]; 
elseif( testGeom == 29 )
    %
    g.dim = 2;
    g.min = [-0.04 - g.dx;0 - g.dx];
    g.max = [1.74 + g.dx; 1.7 + g.dx];
elseif( testGeom == 30 )
    %
    g.dim = 2;
    g.min = [-0.2-g.dx;0 - g.dx];
    g.max = [1.24 + g.dx; 1.24 + g.dx]; 
elseif( testGeom == 31 )
    %
    g.dim = 2;
    g.min = [-0.1-g.dx;0 - g.dx];
    g.max = [1.16 + g.dx; 1.16 + g.dx];
elseif( testGeom >= 32 && testGeom <= 36)
    %
    g.dim = 2;
    h = 0.6;
    g.min = [-h-g.dx;-h - g.dx];
    g.max = [h + g.dx;h + g.dx]; 
elseif( testGeom >= 15 )
    % maskSquarePore.m,
    % maskSquarePore1.m,maskSquarePore2.m,maskSquarePore3.m
    g.dim = 2;
    g.min = [0 - g.dx; 0 - g.dx];
    g.max = [1 + g.dx; 1 + g.dx];
elseif( testGeom == 4) 
    %see maskSphere3D.m 
    g.dim = 3;
    g.min = 0 - g.dx;
    g.max = +1 + g.dx;
elseif( testGeom >= 6 && testGeom <=12)
    % special size for biconic section and throat geometry, gap assumed 0.3 
    %(see maskBiconic2D.m, maskThroat2D.m, maskHybrid1.m
    %     maskThroatBump.m,maskThroatBump1.m, maskThroatHole.m, maskThroat2Dsmall.m)
    g.dim = 2;
    g.min = [-0.5 - g.dx; -0.3 - g.dx];
    g.max = [ 0.5 + g.dx;  0.3 + g.dx];
    
elseif( testGeom == 13 || testGeom == 14)
    % special size for biconic section and throat geometry, gap assumed 0.3 
    %(see maskHybrid2.m, maskHybrid3.m)
    g.dim = 2;
    g.min = [-0.5 - g.dx; -0.4 - g.dx];
    g.max = [ 0.3 + g.dx;  0.4 + g.dx];
else
    g.dim = 2; 
end

%adjust grid for patching, if any
if( testGeom ~= 5 )
    g0 = g;  
    diff = g.max - g.min - [g.dx;g.dx];

    if( patchx || options.addDuct)    
        if( patchx == 3 ) 
            factor = 2;
        else
            factor = 1;
        end
        g.min = g.min - [g.dx;0]
        g.max = g.max + factor*[diff(1);0] +(factor-1)*[g.dx;0];
    end

    if( patchy )
        if( patchy == 3 ) 
            factor = 2;
        else
            factor = 1;
        end

        g.min = g.min - [0;g.dx]
        g.max = g.max + factor*[0;diff(2)]+(factor-1)*[0;g.dx];
    end
end    

% If reading segmented files, basic grid info is taken from seg. data
if(testGeom == 5)
    % reading segmented file name
    [mask,g0] = feval(MaskGeom(testGeom).func,g,options.segFileName); 
    
    % magnify data if desired
    if(options.magFac > 1)
        mask = magnifyData2D(mask,options.magFac);
    end
    
    if(options.addDuct) 
        [m n] = size(mask);
        mask1 = zeros( [m n] );
        for(i=1:n)
                mask1(:,i)=mask(1,i);
        end
        mask = [mask1(1:(m-1),:); mask];
    end
    
    g.min = [1;1]*g.dx;
    [m, n] = size(mask);
    g.max = [m; n]*g.dx;
end



% Create grid
if(periodic)
  if(testGeom ~= 6 ) g.max = (1 - g.dx); end
  g.bdry = @addGhostPeriodic;
else
  g.bdry = @addGhostExtrapolate;
end
g = processGrid(g)


if(patchx || patchy || options.addDuct)
    g0 = processGrid(g0) 
end

% create appropriate mask
if(testGeom ~= 5)
    if (patchx || patchy )
        [mask,g0] = feval(MaskGeom(testGeom).func,g0);
        size(mask)
        [mask1,g0] = feval(MaskGeom(testGeom1).func,g0);
        size(mask1)
        if (patchx > 1 || patchy > 1 ) 
          [mask2,g0] = feval(MaskGeom(testGeom2).func,g0);
          size(mask2)
        end
        
        % patch geometries as specified
        if (patchx == 0)
            mask3 = mask;
            mask4 = mask1;
            if( patchy > 1 )
                mask5 = mask2;
            end
        elseif (patchx == 1)
            mask3 = [mask; mask1(1:(g0.N(1)-1),:)];
            mask4 = [mask1; mask(1:(g0.N(1)-1),:)];
            if (patchy > 1 ) 
              mask5 = [mask2; mask(1:(g0.N(1)-1),:)];
            end
        elseif( patchx == 2 ) 
            mask3 = [mask; mask2(1:(g0.N(1)-1),:)];
            mask4 = [mask2; mask(1:(g0.N(1)-1),:)];
            mask5 = [mask1; mask(1:(g0.N(1)-1),:)];
        elseif( patchx == 3 ) 
            mask3 = [mask;  mask1(1:(g0.N(1)-1),:); mask2(1:(g0.N(1)-1),:)];
            mask4 = [mask2; mask(1:(g0.N(1)-1),:);  mask1(1:(g0.N(1)-1),:)];
            mask5 = [mask1; mask2(1:(g0.N(1)-1),:); mask(1:(g0.N(1)-1),:)];
        end
        
        if(patchy == 0) 
           mask = mask3;
        elseif(patchy == 1)
           mask = [mask3, mask4(:,1:(g0.N(2)-1))];
        elseif(patchy == 2)
           mask = [mask3, mask5(:,1:(g0.N(2)-1))];
        elseif(patchy == 3)
           mask =  [mask3,mask4(:,1:(g0.N(2)-1)),mask5(:,1:(g0.N(2)-1))];
        end
        
        %special case with 4 specified geometries
        if (patchx == 4)
          [mask3,g0] = feval(MaskGeom(testGeom3).func,g0);
          size(mask3)
          
          mask4 = [mask; mask1(1:(g0.N(1)-1),:)];
          mask5 = [mask2; mask3(1:(g0.N(1)-1),:)];
          mask = [mask4, mask5(:,1:(g0.N(2)-1))];
        end  
        size(mask)
        
    elseif(options.addDuct)
        [mask,g0] = feval(MaskGeom(testGeom).func,g0);
        [m n] = size(mask);
        mask1 = zeros( [m n] );
        for(i=1:n)
                mask1(:,i)=mask(1,i);
        end
        mask = [mask1(1:(g0.N(1)-1),:); mask]
        size(mask)
        
    else
        [mask,g] = feval(MaskGeom(testGeom).func,g);
    end
end

if(options.doFlip)
    % flip geometry upside down with respect to the IC
    if(options.ICplane == 'x')
       mask = flipud(mask);
    elseif(options.ICplane == 'y')
       mask = fliplr(mask);
    else
        error('Fix flipping in initalizeDataGrid.m');
    end
end


if(options.doSeal)
   % seal sides of the image according to the IC
    if(options.ICplane == 'x' || options.ICplane == 'c')
       normal = [0.0;-1.0];
       point  = [g.min(1);g.min(2) + 0.5*g.dx(2)];
       normal1 = [0.0;1.0];
       point1  = [g.min(1);g.max(2) - 0.5*g.dx(2)];
       mask1 = max( shapeHyperplane(g,normal,point),shapeHyperplane(g,normal1,point1)); 
       
       mask = max(mask,mask1);
    elseif(options.ICplane == 'y')
       normal = [-1.0;0.0];
       point  = [g.min(1) + 0.5*g.dx(1);g.min(2)];
       normal1 = [1.0;0.0];
       point1  = [g.max(1) - 0.5*g.dx(1);g.min(2)];
       mask1 = max( shapeHyperplane(g,normal,point),shapeHyperplane(g,normal1,point1)); 
       
       mask = max(mask,mask1);
    else
        error('Fix sealing in initalizeDataGrid.m');
    end
end


if(testGeom == 5 || options.doSeal)
    % can't hurt to reinitialize in the case mask is read from a segmnted file
    % also, if any sealing is done, the mask is not a signed distance
    % function any more so it should be reinitialized
    mask = signedDistanceIterative(g,mask,'veryHigh',5*max(g.dx),1e-3);
end

%----------------------------------------------------------------------
% Create initial conditions
%pos = min(g.N(1)/2,g.N(2)/2);
pos = options.IC_pos;
if ( ICplane == 'c' )
    data = circularICx(g,mask);
elseif ( ICplane == 'b' )
    data = bubbleIC(g);
elseif ( ICplane == 'r' )
    data = rectangleIC(g);
else
    data = createIC(g,ICplane,pos);
end

% Need to ensure that the initial conditions satisfy the mask.
data = max(data, mask);
data0 = data;
