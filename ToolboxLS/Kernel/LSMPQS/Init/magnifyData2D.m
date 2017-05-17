function new_data = magnifyData2D(data, mag_factor)

%function magnifyData(data, mag_factor)
% Magnifies data, i.e. replaces each element data(i,j) by 'data(i,j) * ones(mag_factor)'

[m n] = size(data);
m_new = mag_factor * m;
n_new = mag_factor * n;
new_data = zeros( m_new,n_new);

for(i=1:m_new)
  for(j=1:n_new)
      i0 = ceil(i/mag_factor);
      j0 = ceil(j/mag_factor);
      
      new_data(i,j) =  data(i0,j0);
  end
end


