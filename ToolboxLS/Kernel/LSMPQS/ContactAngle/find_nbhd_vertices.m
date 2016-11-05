function  Vnbhd = find_nbhd_vertices(V,F)

% Vnbhd = find_nbhd_vertices(V,F)
% V vertices
% F is a triangulated surfaces (F refers to indices from V)
% Vnbhd has the same number of rows(i.e. vertices) as V
% If c = Vnbhd(v,1) >= 0 then v is boundary vertex and has c neighbors
% whose indices are Vnbhd(v,2:c+1); c is possibly 0 - happens when
% F has some isolated triangles
% Note that "boundary" refers to the essentially 1-dimensional boundary of the
% 2-dimensional triangulated surface F

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
[E I J] = unique(E,'rows');

E_freq = zeros(size(I));

j_sz = size(J); j_sz = j_sz(1);
for(i=1:j_sz)
    E_freq( J(i) ) = E_freq( J(i) ) + 1;
end

Ebdry = E(find(E_freq == 1 ),:);

ebd_sz = size(Ebdry); ebd_sz = ebd_sz(1);
Vind_bdry = reshape(Ebdry,[2*ebd_sz,1]);
Vind_bdry = unique(Vind_bdry);

%Vbdry = V(Vind_bdry,:); %coords for vertices on the boundary of F

%find edges NOT on the boundary of the triangulated surface
Eint = setdiff(E,Ebdry,'rows');

% get number of all vertices
v_num = size(V); v_num = v_num(1);
%indicator array for boundary vertices
Vnbhd = -ones([v_num 1]); % default mark is -1
Vnbhd( Vind_bdry ) = 0; %initial mark for boundary vertices is 0

%use Eint to find interior neighbors of Vbdry
e_sz = size(Eint); e_sz = e_sz(1);

for(i=1:e_sz)
   v1 = Eint(i,1); v2 = Eint(i,2);
   if (Vnbhd(v1) >= 0)
       c = Vnbhd(v1,1) + 1; %this will be counter of neighbors
       Vnbhd(v1,1) = c;
       Vnbhd(v1,c+1) = v2;
   elseif (Vnbhd(v2) >= 0)
       c = Vnbhd(v2,1) + 1; %this will be counter of neighbors
       Vnbhd(v2,1) = c;
       Vnbhd(v2,c+1) = v1;
   end
end
