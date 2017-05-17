function Vbdry = find_bdry_vertices(V,F)

% Vbdry = find_bdry_vertices(V,F)
% creates list of boundary vertices from faces list F (assumed triangles) 
% by exploring F's boundary edges

f_sz = size(F); f_sz = f_sz(1);

% E = list of non-oriented edges, E(i,1) and E(i,2) are vertices on the
% edge sorted such that E(i,1) < E(i,2)
count = 0;
for(i=1:f_sz)
    count = count+1;
    E(count,1) = min(F(i,1),F(i,2));
    E(count,2) = max(F(i,1),F(i,2));
    
    count = count+1;
    E(count,1) = min(F(i,2),F(i,3));
    E(count,2) = max(F(i,2),F(i,3));
    
    count = count+1;
    E(count,1) = min(F(i,3),F(i,1));
    E(count,2) = max(F(i,3),F(i,1));
end

% E has some multiple entries - all repeated ones are internal edges
% all that show up only once are boundary edges.
[Euniq I J] = unique(E,'rows');

Euniq_freq = zeros(size(I));

j_sz = size(J); j_sz = j_sz(1);
for(i=1:j_sz)
    Euniq_freq( J(i) ) = Euniq_freq( J(i) ) + 1;
end

Ebdry = Euniq(find(Euniq_freq == 1 ),:);

ebd_sz = size(Ebdry); ebd_sz = ebd_sz(1);
Vind = reshape(Ebdry,[2*ebd_sz,1]);
Vind = unique(Vind);

Vbdry = V(Vind,:);
