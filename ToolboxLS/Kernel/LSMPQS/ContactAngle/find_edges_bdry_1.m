function Ebdry = find_edges_bdry_1(F)

% Ebdry = find_edges_bdry(F)
% F list of faces(triangles)
% Ebdry = list of non-oriented boundary edges, Ebry(i,1) and Ebdry(i,2) are vertices on the
% edge sorted such that Ebdry(i,1) < Ebdry(i,2)
% Ebdry(i,3) records the face edge i belongs to
% Note that "boundary" refers to the essentially 1-dimensional boundary of the
% 2-dimensional triangulated surface F

E = find_edges_all_1(F);

% E has some multiple entries - all repeated ones are internal edges
% all that show up only once are boundary edges.

E1 = E(:,1:2); %forget third column for now
[E1 I J] = unique(E1,'rows');

E_freq = zeros(size(I));

j_sz = size(J); j_sz = j_sz(1);
for(i=1:j_sz)
    E_freq( J(i) ) = E_freq( J(i) ) + 1;
end

E2 = E(I,:);  % PUT THE THIRD COLUMN BACK
Ebdry = E2(find(E_freq == 1 ),:);

