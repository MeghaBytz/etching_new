function [angle1 angle2 angle3 angle4] = run_test

cd /ices/masha/sept00_sub256/swr/fld_seg/
tic
angle1 = explore_contact_angle_ubc('w_swr_ubc.gz','o_swr_ubc.gz',0,256,256,256);
toc

cd /ices/masha/sept00_sub256/sor/fld_seg/
tic
angle2 = explore_contact_angle_ubc('w_sor_ubc.gz','o_sor_ubc.gz',0,256,256,256);
toc


cd /ices/masha/sept00_sub256/gelswr/fld_seg/
tic
angle3 = explore_contact_angle_ubc('w_gelswr_ubc.gz','o_gelswr_ubc.gz',0,256,256,256);
toc


cd /ices/masha/sept00_sub256/gelsor/fld_seg/
tic
angle4 = explore_contact_angle_ubc('w_gelsor_ubc.gz','o_gelsor_ubc.gz',0,256,256,256);
toc