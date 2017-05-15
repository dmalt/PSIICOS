clear all;
protocol_path = '/home/dmalt/fif_matlab/Brainstorm_db/motorprobe';
subj_ID = 'AB';
suffix = '_control_RH_raw_tsss_ica';

condition = [subj_ID, suffix]; 
freq_band = [10, 20];
time_range_pre = [-0.7, -0.2];
time_range_post = [0, 0.5];

GainSVDTh = 0.01;
isLR = true;
lambda = 1;
pwr_rnk = 350;
Upwr = [];
threshold = 100;

[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subj_ID, protocol_path);

HM = ups.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);

trials_pre = ups.LoadTrials(subj_ID, condition,...
                                freq_band, time_range_pre,...
                                HM, protocol_path);

trials_post = ups.LoadTrials(subj_ID, condition,...
                                 freq_band, time_range_post,...
                                 HM, protocol_path);

tr_pre = trials_pre.data;%(:, :, 51:100);% + 50 * rand(size(trials.data));
tr_post = trials_post.data;%(:, :, 51:100);% + 50 * rand(size(trials.data));

CT_pre = ups.CrossSpectralTimeseries(tr_pre, true);
CT_post = ups.CrossSpectralTimeseries(tr_post, true);

CT_resh_pre = reshape(mean(CT_pre, 2), sqrt(size(CT_pre, 1)), sqrt(size(CT_pre, 1)));
CT_resh_post = reshape(mean(CT_post, 2), sqrt(size(CT_post, 1)), sqrt(size(CT_post, 1)));

[~, Ps_pre] = ups.DICS(CT_resh_pre, HM.gain, lambda, false);
[~, Ps_post] = ups.DICS(CT_resh_post, HM.gain, lambda, false);


figure;
h = plot_brain_cmap(CtxHHR, Ctx, [], Ps, zeros(size(Ps)), 0.2);

CT_proj_from_vc = ps.ProjectAwayFromPowerFixedOr(CT_resh(:), HM.gain, pwr_rnk);
[CS, IND] = ps.PSIICOS_ScanFast(HM.gain, (CT_proj_from_vc));

% CT_rand = rand(size(CT_resh)) *  + j * rand(size(CT_resh));
% [Cs, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_rand), HM.gain, pwr_rnk, Upwr);
% [CS, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_resh), HM.gain, pwr_rnk, Upwr);

con_inds = ups.threshold_connections(CS, threshold, IND);

con = ups.Connections(subj_ID, con_inds, freq_band, time_range, CT, condition, HM, CtxHHR);
% con_c = con.Clusterize(10,0.02);
figure;
con.Plot(0.2);
