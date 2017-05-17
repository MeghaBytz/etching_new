function fun = correct_fld_intfc_bdry(fun,Vnbhd)
% fun = correct_fld_intfc_bdry(fun,Vnbhd)
% Correct interpolated function 'fun' values at fld intfc boundary points.
% Vnbhd is supposed to be returned by function fld_nbhd_vertices()

Vbdry_ind = find( Vnbhd(:,1) > 0 );

v_sz = size(Vbdry_ind); v_sz = v_sz(1);

for(i=1:v_sz)
    j = Vbdry_ind(i);
    c = Vnbhd(j,1); % number of neighbors for the boundary node
    new_val = 0;
    for(k=2:(c+1))
        new_val = new_val + fun( Vnbhd(j,k) ); 
    end
    % replace value at j'th node with average of neighbor's values
    fun(j) = new_val / c; 
end
