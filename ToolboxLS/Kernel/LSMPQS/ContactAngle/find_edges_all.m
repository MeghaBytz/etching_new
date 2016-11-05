function E = find_edges_all(F)

% E = find_edges_all(F)
% F list of faces(triangles)
% E = list of non-oriented edges, E(i,1) and E(i,2) are vertices on the
% edge sorted such that E(i,1) < E(i,2)
% the list HAS REPETITIONS

f_sz = size(F); f_sz = f_sz(1);

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
