clear all;
%  ----------------------------------- Simulations ---------------- %
% subjID = 'test';
% freqBand = [];
% t_range = [];
% [HM, CT, trials, Ctx, XYZGenOut] = ups.SimulateData(pi/2, 100, 0.01, 0.5, 0, false);
% ----------------------------------------------------------------- %


% ----------------------- Real data -------------------------------- %
% finger
% protocol_path = '/home/dmalt/fif_matlab/Brainstorm_db/motorprobe';
% subj_ID = 'AK';
% suffix = '_control_RH_raw_tsss_ica';

% condition = [subj_ID, suffix]; 
% freq_band = [10, 20];
% time_range = [-0.9, -0.5];
% GainSVDTh = 0.01;
% isLR = true;
% sFreq  = 1000;

% odd-ball
protocol_path = '~/PSIICOS_osadtchii';
subj_ID = '0109_zvma';
condition = '2';
freq_band = [16;25];
time_range = [0;0.7];
sFreq = 500;
isLR = true;
GainSVDTh = 0.01;

trials = ups.LoadTrials(subj_ID, condition,...
                        freq_band, time_range, sFreq,...
                        GainSVDTh, protocol_path);
trials = trials.data;

HM = ups.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);
[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subj_ID, protocol_path);
% ---------------------------------------------------------------- %

solution = [];
Upwr = [];
rel_threshold = 10;
n_resamples = 100;
pwr_rnk = 500;

CT_true = ups.CrossSpectralTimeseries(trials);

[INDtrue, Cp, Upwr, Cs_true, IND] = ps.T_PSIICOS(CT_true, HM.gain, rel_threshold, pwr_rnk, 0, Upwr);
max_cs_true = max(Cs_true);

max_cs_surr = zeros(n_resamples, 1);

for i = 1:n_resamples
    % [spoilt_trials, solution] = ups.ShufflePhases(trials, i, size(HM.gain,1), solution);
    [spoilt_trials, solution] = ups.ShufflePhases(trials, i, 30, solution);
    CT = ups.CrossSpectralTimeseries(spoilt_trials);
    [INDsurr{i}, Cp, Upwr, Cs, IND] = ps.T_PSIICOS(CT, HM.gain, rel_threshold, pwr_rnk, 0, Upwr);
    max_cs_surr(i) = max(Cs)
end

p_val = sum(max_cs_surr > max_cs_true) / n_resamples;
disp(p_val);

con_ps_surr = ups.Connections(subj_ID, INDsurr,  freq_band, time_range, CT, 'test', HM, CtxHHR);
con_ps_av = con_ps_surr.Average();
con_ps_true = ups.Connections(subj_ID, INDtrue,  freq_band, time_range, CT, 'test', HM, CtxHHR);



CT_resh = reshape(mean(CT,2), sqrt(size(CT,1)), sqrt(size(CT,1)));
CT_resh_true = reshape(mean(CT_true,2), sqrt(size(CT_true,1)), sqrt(size(CT_true,1)));
[A, Ps] = ups.DICS(CT_resh, HM.gain, 10);
[A_true, Ps_true] = ups.DICS(CT_resh_true, HM.gain, 10);

figure;
h = plot_brain_cmap(CtxHHR, Ctx, [], Ps, zeros(size(Ps)), 0.2);
figure;
h = plot_brain_cmap(CtxHHR, Ctx, [], Ps_true, zeros(size(Ps)), 0.2);
% [SPC, TPR, PPV] = ups.GenerateROC(Cs, 0.02, HM.GridLoc, IND, 200, XYZGenOut, 1);
