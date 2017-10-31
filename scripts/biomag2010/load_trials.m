% Load trials for mental rotation dataset
% --------------------------------------- %
% AUTHOR: dmalt
% DATE: Tue Oct 31 14:26:03 MSK 2017
% --------------------------------------- %


% -------- setup  constants ------- %

% paths
home = getenv('HOME');
data_path = [home, '/Data/mentrot/MentalRotationDeLange/preprocessed'];
subj_ID = 'biomag2010';
protocol_path = [home, '/bst/brainstorm_db/mentrot/'];
condition = 'raw';

% head model
isLR = true;
GainSVDTh = 0.01;

% trials
time_range = [0,1];
freq_band = [8,12];

% --------------------------------- %

% load trials
data_files = dir([data_path, '/dataset*']);
data_fnames =  {data_files.name};

i_file = 1;
load_fname = [data_path, '/', data_fnames{i_file}];
data = load(load_fname);
data = data.data;
trials = data.trial;
s_freq = data.fsample;
times = data.time;
n_times = length(find((times{1} >= time_range(1)) & (times{1} < time_range(2))));


% load head model
HM = ups.bst.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);

% filter trials
n_ch = size(HM.gain, 1);
n_tr = length(trials);
tr_filt = zeros(n_ch, n_times, n_tr);
[b,a] = fir1(128, freq_band / (s_freq / 2), 'bandpass');
for i_tr = 1:length(trials)
    tmp = filtfilt(b, a, (HM.UP * trials{i_tr})')'; % filter and reduce dim
    samp_range =  find((times{i_tr} >= time_range(1)) & (times{i_tr} < time_range(2)));
    tr_filt(:,:, i_tr) = 1e12 * tmp(:, samp_range); % crop
end


% compute cross-spectrum timeseries
CT = ups.conn.CrossSpectralTimeseries(tr_filt, true);
CT_resh = reshape(mean(CT,2), sqrt(size(CT,1)), sqrt(size(CT,1)));
lambda = 0.2;
[A, Ps] = ups.conn.DICS(CT_resh, HM.gain, lambda);
[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subj_ID, protocol_path);
figure;
h = plot_brain_cmap(CtxHHR, Ctx, [], Ps, zeros(size(Ps)), 0.2);

