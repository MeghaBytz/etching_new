function angle = calculate_angle_vec(N1,N2,do_vis)

%angle = calculate_angle_vec(N1,N2,do_vis)
% N1,N2 assumed normalized vectors

    prod = N1 .* N2;
    
    cos_angle =( prod(:,1) + prod(:,2) + prod(:,3) );
    
    % results > 1.0 or < -1.0 are rounding off errors, correct for those
    cos_angle( find( cos_angle > 1.0) ) = 1.0;
    cos_angle( find( cos_angle < -1.0) ) = -1.0;

    angle = acos(cos_angle);
    angle = angle * (180 / pi);

    if( do_vis ) figure, hist(angle,50);  end