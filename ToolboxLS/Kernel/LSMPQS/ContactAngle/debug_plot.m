function debug_plot(infl_name1,infl_name2)

% debug_plot(infl_name1,infl_name2)
% visualizes triple points found in 'explore_contact_angle' as blue dots
% over the fluid1-fluid2 interface surface plot in gray

[V angle_intrp] = explore_contact_angle(infl_name1,infl_name2);

figure, h_fld = patch('Faces',F_FLD,'Vertices',V1);
set(h_fld,'FaceColor','y','EdgeColor', 'none');
daspect([1 1 1]); camlight; lighting phong;
view(3), hold on

% Note - grid minimum g.min in explore_curv_fld() is 0 and in
% explore_contact_angle() is 1, so we need to adjust in order to plot
% vertices on the same plot
V = V - ones(size(V));

cut_off1 = 40; cut_off2 = 90;
i = find( angle_intrp < cut_off1);
j = find( (angle_intrp > cut_off1 )  & (angle_intrp < cut_off2 ) );
k = find(angle_intrp > cut_off2);

plot3(V(i,1),V(i,2),V(i,3),'b.')
plot3(V(j,1),V(j,2),V(j,3),'g.')
plot3(V(k,1),V(k,2),V(k,3),'r.')