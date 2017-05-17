function Ebdry = find_edges_bdry(F)

% Ebdry = find_edges_bdry(F)
% F list of faces(triangles)
% Ebdry = list of non-oriented boundary edges, Ebry(i,1) and Ebdry(i,2) are vertices on the
% edge sorted such that Ebdry(i,1) < Ebdry(i,2)
% Note that "boundary" refers to the essentially 1-dimensional boundary of the
% 2-dimensional triangulated surface F

E = find_edges_all(F);

% E has some multiple entries - all repeated ones are internal edges
% all that show up only once are boundary edges.
[E I J] = unique(E,'rows');

E_freq = zeros(size(I));

j_sz = size(J); j_sz = j_sz(1);
for(i=1:j_sz)
    E_freq( J(i) ) = E_freq( J(i) ) + 1;
end

Ebdry = E(find(E_freq == 1 ),:);

