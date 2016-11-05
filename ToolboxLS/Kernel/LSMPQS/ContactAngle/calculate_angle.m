function angle = calculate_angle(dx_mask,dy_mask,dz_mask,dx_data,dy_data,dz_data,do_vis)

    prod = dx_mask .* dx_mask + dy_mask .* dy_mask + dz_mask .* dz_mask;
    prod1 = dx_data .* dx_data + dy_data .* dy_data + dz_data .* dz_data;

    prod1 = prod .* prod1;  prod1 = sqrt(prod1);
    prod1( find( prod1 == 0) ) = 1; % avoid division by zero

    prod3 = dx_mask .* dx_data + dy_mask .* dy_data + dz_mask .* dz_data;
    cos_angle = prod3 ./ prod1; 
    % results > 1.0 or < -1.0 are rounding off errors, correct for those
    cos_angle( find( cos_angle > 1.0) ) = 1.0;
    cos_angle( find( cos_angle < -1.0) ) = -1.0;

    angle = acos(cos_angle);
    angle = angle * (180 / pi);

    if( do_vis ) figure, hist(angle,50);  end

    num_of_triple_points = size(angle)
    median_angle  = median(angle)
