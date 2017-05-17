function  [ADJ_FLD ADJ_GRAIN ] = find_adjacent_faces(F_FLD,F_GRAIN)

% [ADJ_FLD ADJ_GRAIN ] = find_adjacent_faces(F_FLD,F_GRAIN)
% F_FLD fluid1-fluid2 faces(triangles)
% F_GRAIN fluid1-grain faces(triangles)
% ADJ_FLD and ADJ_GRAIN are indices of adjacent F_FLD and F_GRAIN
% triangles (i.e. the ones that share an edge)

Ebdry = find_edges_bdry_1(F_FLD);
% For each bdry edge, the face it belongs to is recorded in Ebdry(:,3)
ADJ_FLD = Ebdry(:,3);

Egrain = find_edges_all_1(F_GRAIN);

e_sz = size(Ebdry); e_sz = e_sz(1);

for(i=1:e_sz)
   v1 = Ebdry(i,1); v2 = Ebdry(i,2); % note that  v1 < v2 (always)
   % looking for entries in Eg1 such that Eg1(j,1) == v1 and Eg1(j,2) == v2
   A = ismember(Egrain,v1); A = A(:,1);
   B = ismember(Egrain,v2); B = B(:,2);
   J = find( A .* B );
   j_sz = size(J); j_sz = j_sz(1);
   if( j_sz ) 
       for(k=1:j_sz)
           ADJ_GRAIN(i,k) = Egrain(J(k),3);
       end
   else
       ADJ_GRAIN(i,1) = 0; % adjacent grain traingle is not found
   end
end

valid = find( ADJ_GRAIN(:,1) > 0 );
ADJ_GRAIN = ADJ_GRAIN(valid);
%ADJ_GRAIN = ADJ_GRAIN(valid,:); %uncoment if more than one adjacent grain
%edge makes sense in the context
ADJ_FLD   = ADJ_FLD(valid);