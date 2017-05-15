clear all;
protocol_path = '/home/dmalt/fif_matlab/Brainstorm_db/motorprobe';
subj_ID = 'AB';
suffix = '_control_RH_raw_tsss_ica';

condition = [subj_ID, suffix]; 
freq_band = [19, 25];
time_range_pre = [-0.7, -0.2];
time_range_post = [0, 0.5];

GainSVDTh = 0.01;
isLR = true;
lambda = 1;
SL_rnk = 350;
sig_rnk = 0;
cp_part = 'real';
is_fast = false;
Upwr = [];
threshold = 100;
seed_ind = 1;

[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subj_ID, protocol_path);

HM = ups.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);
[tr_pre, CT_pre, Ps_pre, CS_pre] = prepare_cond(subj_ID, condition, protocol_path,...
                                        HM, isLR, GainSVDTh, freq_band,...
                                        time_range_pre, lambda, seed_ind,...
                                        SL_rnk, sig_rnk, Upwr, cp_part, is_fast);

[tr_post, CT_post, Ps_post, CS_post] = prepare_cond(subj_ID, condition, protocol_path,...
                                        HM, isLR, GainSVDTh, freq_band,...
                                        time_range_post, lambda, seed_ind,...
                                        SL_rnk, sig_rnk, Upwr, cp_part, is_fast);
figure;
h = plot_brain_cmap(CtxHHR, Ctx, [], Ps_pre, zeros(size(Ps_pre)), 0.2);

figure;
h = plot_brain_cmap_hemisplit(CtxHHR, Ctx, [], CS_pre.data, zeros(size(CS_pre.data)), 0, [0,0,0])
% CT_proj_from_vc = ps.ProjectAwayFromPowerFixedOr(CT_resh(:), HM.gain, pwr_rnk);
% [CS, IND] = ps.PSIICOS_ScanFast(HM.gain, (CT_proj_from_vc));

% CT_rand = rand(size(CT_resh)) *  + j * rand(size(CT_resh));
% [Cs, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_rand), HM.gain, pwr_rnk, Upwr);
% [CS, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_resh), HM.gain, pwr_rnk, Upwr);

% con_inds = ups.threshold_connections(CS, threshold, IND);

% con = ups.Connections(subj_ID, con_inds, freq_band, time_range, CT, condition, HM, CtxHHR);
% con_c = con.Clusterize(10,0.02);
% figure;
% con.Plot(0.2);
