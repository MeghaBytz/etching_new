function process_blobs_contact_angle(basename,blob_id_list,phase1,phase2)

% process_blobs_contact_angle(basename,blob_id_list,phase1,phase2)
% basename - blob fluid segmented files basename w/o blob id or fluid
% identifier; files assumed gzipped
%           -- basically a basename given to blob segfiles in case7.14,
%           3DMA-Rock
% blob_id_list - files storing id's of blobs to process
% phase1,2 - phase identifiers e.g. OIL,H2O
% if processing oil blobs put OIL as first identifier and vice versa

fid = fopen(blob_id_list,'r');
[blob_id blob_count] = fscanf(fid,'%d\n');

do_vis = 0;
total_sum = 0; total_num = 0;
for(i=1:blob_count)
   filename1 = sprintf('%s_%d_%s_sgn_dist.gz',basename,blob_id(i),phase1);
   filename2 = sprintf('%s_%d_%s_sgn_dist.gz',basename,blob_id(i),phase2);
   
   [angle median_angle1 mean_angle1] = explore_contact_angle(filename1,filename2,do_vis);
   median_angle(i) = median_angle1;
   mean_angle(i) = mean_angle1;
   
   total_sum = total_sum + sum(angle);
   tmp_sz = size(angle); tmp_sz = tmp_sz(1);
   total_num = total_num + tmp_sz;
end

save median_angle median_angle
save mean_angle mean_angle
total_mean = total_sum/total_num;
save total_mean total_mean