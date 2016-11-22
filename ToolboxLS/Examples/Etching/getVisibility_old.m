function vis = getVisibility_old(g, source, data, mask, options)

intfcPoints = isNearInterface(data);
w = zeros(size(data,1),1);
h = zeros(size(data,2),1);
intfcPoints(mask >= 0) = 0;
vis = -1*ones(size(data));

if (options.doMask)
    source_angle = atan(abs((g.xs{1,1}-source(1,1))./(g.xs{2,1} - source(1,2))));
    
    % locate the bottom of the mask, and the hole dimensions.
    intfcPoints_mask = isNearInterface(mask);
    
    a = find(intfcPoints_mask(1,:));
    mask_bot = a(2); mask_top = a(end-1);
    
    b = find(intfcPoints_mask(:,mask_bot)==0);
    left_bot_corner = [g.min(1) + b(1)*g.dx(1), g.min(2) + mask_bot*g.dx(2)]; 
    right_bot_corner = [g.min(1) + b(end)*g.dx(1), g.min(2) + mask_bot*g.dx(2)];

    c = find(intfcPoints_mask(:,mask_top) == 0);
    left_top_corner = [g.min(1) + c(1)*g.dx(1), g.min(2) + mask_top*g.dx(2)]; 
    right_top_corner = [g.min(1) + c(end)*g.dx(1), g.min(2) + mask_top*g.dx(2)];
    
    angle_left = atan((source(1,1) - left_bot_corner(1,1))/(source(1,2) - left_bot_corner(1,2)));
    angle_right = atan(right_bot_corner(1,1) - (source(1,1))/(source(1,2) - right_bot_corner(1,2)));
    
    vis((g.xs{1,1} < source(1,1)) & (source_angle > angle_left)) = 0;
    vis((g.xs{1,1} > source(1,1)) & (source_angle > angle_right)) = 0;
    vis(data < 0) = 0;
    
    for i = 1:201
        for j = 1:201
            if(~vis(i,j))
                for k = j:201
                    for l = 1:201
                        if((vis(l,k)==0) && (source_angle(l,k) < source_angle(i,j)))
                            vis(i,j) = 0;
                        end
                    end
                end
            end
        end
    end
    
    

end