% Load forward operator and timeseries
% fname = 'K0021_data.mat';
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
freqBand = [18, 22];
sFreq = 500;                             % can I figure this out from the data?
[b,a] = fir1(128, freqBand / (sFreq / 2), 'bandpass');    % define filter
TS_pca_filt = filtfilt(b, a, TS_pca');
TS_pca_filt_hilb = hilbert(TS_pca_filt)';
% ------------------------------------------------------ %

window_size = 500;
window_step = 100;
% - Slide along the recording with a window and compute time-averaged cross-spectrum - %
if ~exist('CP_time_av.mat', 'file')
    pwr_rnk = 350;
    window = 1:window_size;
    CP_time_av = zeros(n_sensors_reduced ^ 2, fix(n_samples / window_step) + 1);

    for iTime = 1:fix(n_samples / window_step)
        for jTime = window
            tmp = TS_pca_filt_hilb(:, jTime) * TS_pca_filt_hilb(:, jTime)'; 
            CP_time_av(:, iTime) = CP_time_av(:, iTime) + tmp(:);  
        end
        
        CP_time_av(:, iTime) = CP_time_av(:, iTime) / window_size; 
        CP_time_av(:, iTime) = ProjectAwayFromPowerComplete(CP_time_av(:,iTime), G_tang_pca, pwr_rnk);
        if window(end) + window_step < n_samples
            window = window + window_step;
        else 
            break
        end
    end
    save('CP_time_av.mat', 'CP_time_av');
else
    load('CP_time_av.mat');
end
% ------------------------------------------------------------------------------------ %
if ~exist('pairs.mat', 'file')
    for iTime = 1:size(CP_time_av, 2)
        [Cs, IND] = MUSIC_ScanFast(G_tang_pca, CP_time_av(:,iTime));
        pairs{iTime} = get_ij_above_threshold(Cs, 0.8, IND);
    end
    save('pairs.mat', 'pairs');
else
    load('pairs.mat');
end

% iWindow = 110
if ~exist('con_av.mat','file')
    for iWindow = 1:length(pairs)
        iWindow
        con = Bundles(pairs{iWindow}(1:20:end,:), HM, CtxHR );

        con_clust = con.Clusterize(20, 0.02);
        con_av{iWindow} = con_clust.Average();
        % con_av.Plot()
    end
    % save('con_av.mat', 'con_av', '-v7.3');
else
    load('con_av.mat');
end


