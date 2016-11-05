function angle = calculate_face_angles(V,F_FLD_INTFC,ADJ_FLD,F_GRAIN,ADJ_GRAIN)

% angle = calculate_face_angles(V,F_FLD_INTFC,ADJ_FLD,F_GRAIN,ADJ_GRAIN
% V - vertices used by faces in F_FLD_INTFC and F_GRAIN
% ADJ_FLD and ADJ_GRAIN are arrays of faces IDs that correspond to faces
% that share an edge
% note that faces are actually triangles
% angle(i) is an angle between triangle ADJ_FLD(i) and ADJ_GRAIN(i)

num = size(ADJ_FLD); num = num(1)

for(f=1:num)
   ff =  ADJ_FLD(f); 
   %find triangle vertices for face ff in F_FLD_INTFC
   v(1) = F_FLD_INTFC(ff,1);  v(2) = F_FLD_INTFC(ff,2); v(3) = F_FLD_INTFC(ff,3); 
   
   fg = ADJ_GRAIN(f);
   if( fg > 0 )
       %find triangle vertices for face fg in F_GRAIN
       w(1) = F_GRAIN(fg,1);   w(2) = F_GRAIN(fg,2);    w(3) = F_GRAIN(fg,3); 

       % no need to reorder since I can't trust that faces from 'isosurface' are output in a
       % consistently oriented way
       %[c i j] = intersect(v,w); 
       %i(3) = setdiff([1 2 3],i);
       %j(3) = setdiff([1 2 3],j);
       %v = v(i);  w = w(j);  % reorder vertices

       p = V(v,:);   r = V(w,:);

       %calculate normal to the fld intfc triangle plane
       n1 = cross(p(2,:)-p(1,:),p(3,:)-p(2,:));
       norm1 = norm(n1);
       if (norm1 < 1000*eps) disp('Bummer - calculate_face_angles()'); end

       %calculate normal to the grain intfc triangle plane
       n2 = cross(r(2,:)-r(1,:),r(3,:)-r(2,:));
       norm2 = norm(n2);
       if (norm2 < 1000*eps) disp('Bummer - calculate_face_angles()'); end

       %calculate angle btw normals
       cos_angle(f) = dot(n1,n2) / (norm1*norm2);
   else
       cos_angle(f) = 1000.0 %mark invalid entries
   end
   
end

%get rid of invalid entries
cos_angle = cos_angle( find (cos_angle ~= 1000.0 ));
% get rid of rounding errors
cos_angle( find( cos_angle > 1.0) ) = 1.0;
cos_angle( find( cos_angle < -1.0) ) = -1.0;

angle = acos(cos_angle);
angle = angle * (180 / pi);
% do not look for minimum value any more
%angle = min(angle,180-angle);