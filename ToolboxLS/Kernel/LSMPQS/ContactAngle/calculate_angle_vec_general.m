function angle = calculate_angle_vec_general(N1,N2,do_vis)

%angle = calculate_angle_vec(N1,N2,do_vis)
% N1,N2 assumed NOT normalized vectors

    size(N1)
    size(N2)
    prod = N1 .* N2;
    
    norm1 = N1 .* N1; norm1 = sqrt(norm1(:,1) + norm1(:,2) + norm1(:,3));
    norm2 = N2 .* N2; norm2 = sqrt(norm2(:,1) + norm2(:,2) + norm2(:,3));
    norm = norm1 .* norm2;
    norm ( find(norm == 0 )) = 1;
    
    cos_angle =( prod(:,1) + prod(:,2) + prod(:,3) ) ./ norm;
    
    % results > 1.0 or < -1.0 are rounding off errors, correct for those
    cos_angle( find( cos_angle > 1.0) ) = 1.0;
    cos_angle( find( cos_angle < -1.0) ) = -1.0;

    angle = acos(cos_angle);
    angle = angle * (180 / pi);

    if( do_vis ) figure, hist(angle,50);  end

    num_of_vertices = size(angle)
    median_angle  = median(angle)