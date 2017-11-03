function tr_filt = load_trials(data_path, time_range, freq_band, HM)
% Load trials for mental rotation dataset
% --------------------------------------- %
% AUTHOR: dmalt
% DATE: Tue Oct 31 14:26:03 MSK 2017
% --------------------------------------- %

% load trials
    data_files = dir([data_path, '/dataset*']);
    data_fnames =  {data_files.name};

    times = [];
    trials = [];
    for i_file =1:5
        load_fname = [data_path, '/', data_fnames{i_file}];
        data = load(load_fname);
        data = data.data;
        trials = [trials, data.trial];
        times = [times, data.time];
    end
    s_freq = data.fsample;
    n_times = length(find((times{1} >= time_range(1)) & (times{1} < time_range(2))));


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
end
