ch_used = 36:186;

data_path = '/home/dmalt/Documents/MATLAB/bst/brainstorm_db/mentrot/data/biomag2010/@rawdataset02';
hm_hr_path = [data_path, '/headmodel_surf_os_meg.mat'];
hm_lr_path = [data_path, '/headmodel_surf_os_meg_02.mat'];

hm_hr = load(hm_hr_path);
hm_lr = load(hm_lr_path);


G_hr = hm_hr.Gain(36:186,:);
G_lr = hm_lr.Gain(36:186,:);

G2_hr = ups.ReduceToTangentSpace(G_hr);
G2_lr = ups.ReduceToTangentSpace(G_lr);

GridLoc_hr = hm_hr.GridLoc;
GridLoc_lr = hm_lr.GridLoc;

save_fname = 'ctf_hm_hr_lr_with_gridloc_full_dim.mat';

save(save_fname, 'G2_hr', 'G2_lr', 'GridLoc_hr', 'GridLoc_lr', '-v7.3');
