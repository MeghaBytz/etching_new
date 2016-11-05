function [ind comp size_comp] = findConnComp(data,dim,pos)

% function [ind comp size_comp] = findConnComp(data,dim,pos)
% data - 2d array
% dim - dimension 1 or 2
% pos - positions (in specified dimension)
% finds connected components the negative 'data' values in column (row)
% 'pos' - column  or row depends on 'dim' choice
% Function returns
% ind - indices of all negative 'data' values
% comp - corresponding connected component (1,2,...)
% size_comp - size_comp(i) is size of the i-th connected component

if(dim == 1)
   [i j] = find(data(pos,:) < 0);
   check_ind = j; %note that each 'i' element is equal to 'pos'
   ind = find(data(pos,:) < 0);
   
else
   [i j] = find(data(:,pos) < 0 );
   check_ind = i; %note that each 'j' element is equal to 'pos'
   ind = find(data(:,pos) < 0);
end

n = numel(ind);
if (n == 0) 
    num_comp = 0;
    comp = 0;
    size_comp = 0;
    return;
end

num_comp = 1; size_comp(num_comp) = 1;
k = 1;
comp(k) = 1;
for(k=2:n)
    if(check_ind(k) == check_ind(k-1)+1 )
        % we're still in the same component
        comp(k) = num_comp;
        size_comp(num_comp) = size_comp(num_comp) + 1;
    else
        % new component
        num_comp = num_comp+1; size_comp(num_comp) = 1;
        comp(k) = num_comp; 
    end
end
