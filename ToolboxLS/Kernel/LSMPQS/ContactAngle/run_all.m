function run_all

cd /ices/masha/sept00_sub256/swr/blobs/
process_blobs_contact_angle('h2oblob','h2oblob_list_interior','H2O','OIL');

cd /ices/masha/sept00_sub256/gelswr/blobs/
process_blobs_contact_angle('h2oblob','h2oblob_list_interior','H2O','OIL');

cd /ices/masha/sept00_sub256/sor/blobs/
process_blobs_contact_angle('oilblob','oilblob_list_interior','OIL','H2O');

cd /ices/masha/sept00_sub256/gelsor/blobs/
process_blobs_contact_angle('oilblob','oilblob_list_interior','OIL','H2O');
