main
% trials
time_range_pre = [-0.5,-0.2];
time_range_post = [0,1];

beta_band = [16,24];
alpha_band = [8,12];
gamma_band = [30,60];

% freq_band = beta_band;
freq_band = [8,30];
% freq_band = alpha_band;
% freq_band = gamma_band;

ltr = memoize(@load_trials);

tr_filt_pre = ltr(data_path, time_range_pre, freq_band, HM);
tr_filt_post = ltr(data_path, time_range_post, freq_band, HM);

% prepare sensors
ch_path = [protocol_path, '/data/biomag2010/@rawdataset02/channel_ctf_acc1.mat'];
meg_loc_3d  = prepare_sensors(ch_path);

% compute cross-spectrum timeseries
CT_pre = ups.conn.CrossSpectralTimeseries(tr_filt_pre, true);
CT_post = ups.conn.CrossSpectralTimeseries(tr_filt_post, true);

CT_resh_pre = reshape(mean(CT_pre,2), sqrt(size(CT_pre,1)), sqrt(size(CT_pre,1)));
CT_resh_post = reshape(mean(CT_post,2), sqrt(size(CT_post,1)), sqrt(size(CT_post,1)));
pwr_pre = HM.UP' * diag(CT_resh_pre);
pwr_post = HM.UP' * diag(CT_resh_post);

delta = pwr_post - pwr_pre;% ./ pwr_pre;
ratio = zeros(size(delta));
ratio(delta > 0) = pwr_post(delta > 0) ./ pwr_pre(delta > 0);
ratio(delta < 0) = - pwr_pre(delta < 0) ./ pwr_post(delta < 0);

ratio(ratio > 3) = 3;
ratio(ratio < -3) = -3;
plot_topo(meg_loc_3d, ratio);
colorbar;
caxis([-2,2]);
