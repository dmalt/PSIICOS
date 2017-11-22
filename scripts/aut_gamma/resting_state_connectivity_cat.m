% Load forward operator and timeseries
fname = 'timeseries_K0021_eo.mat';
load(fname)
n_samples = size(timeseries, 2);
n_sensors = size(timeseries, 1);

% G_tang = ReduceToTangentSpace(leadfield, 'all');
% [G_tang_pca, UP] = ReduceSensorSpace(G_tang, 0.01);

HM = LoadHeadModel('K0021', 'K0021_rest_raw_tsss_mc_trans_ec', '/home/dmalt/fif_matlab/Brainstorm_db/aut_gamma/');
[Ctx, CtxHR] = GetCtx('K0021', '/home/dmalt/fif_matlab/Brainstorm_db/aut_gamma/');
TS_pca = HM.UP * timeseries;
G_tang_pca = HM.gain;
n_sensors_reduced = size(HM.UP,1);


% ---- Bandpass filter + analytic signal computation --- %
freqBand = [16,25];
sFreq = 500;                             % can I figure this out from the data?
[b,a] = fir1(128, freqBand / (sFreq / 2), 'bandpass');    % define filter
TS_pca_filt = filtfilt(b, a, TS_pca');
TS_pca_filt_hilb = hilbert(TS_pca_filt)';
% ------------------------------------------------------ %

% - Slide along the recording with a window and compute time-averaged cross-spectrum - %
pwr_rnk = 750;
window_size = 500;
window_step = 100;
window = 1:window_size;
CP_time_av = zeros(n_sensors_reduced ^ 2, 1);

for jTime = 1:n_samples
	tmp = TS_pca_filt_hilb(:, jTime) * TS_pca_filt_hilb(:, jTime)'; 
	CP_time_av = CP_time_av + tmp(:);  
end

CP_time_av = CP_time_av / n_samples; 
CP_time_av = ProjectAwayFromPowerComplete(CP_time_av, G_tang_pca, pwr_rnk);
% ------------------------------------------------------------------------------------ %
[Cs, IND] = MUSIC_ScanFast(G_tang_pca, real(CP_time_av));
pairs = get_ij_above_threshold(Cs, 0.3, IND);

% iWindow = 110
con = Bundles(pairs(1:20:end,:), HM, CtxHR );
con_clust = con.Clusterize(10, 0.02);
con_av = con_clust.Average();
