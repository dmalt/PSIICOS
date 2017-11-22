clear all;
protocol_path = '/home/dmalt/fif_matlab/Brainstorm_db/motorprobe';
subj_ID = 'AP';
suffix = '_control_RH_raw_tsss_ica';

condition = [subj_ID, suffix]; 
freq_band = [10, 20];
time_range = [-0.7, -0.2];
GainSVDTh = 0.01;
isLR = true;
lambda = 30;
pwr_rnk = 350;
Upwr = [];
threshold = 100;

trials = ups.LoadTrials(subj_ID, condition,...
                        freq_band, time_range,...
                        GainSVDTh, protocol_path);

HM = ups.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);

tr = trials.data;%(:, :, 51:100);% + 50 * rand(size(trials.data));
CT = ups.CrossSpectralTimeseries(tr, true);
CT_resh = reshape(mean(CT,2), sqrt(size(CT,1)), sqrt(size(CT,1)));
[A, Ps] = ups.DICS(CT_resh, HM.gain, lambda);
[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subj_ID, protocol_path);

figure;
h = plot_brain_cmap(CtxHHR, Ctx, [], Ps, zeros(size(Ps)), 0.2);

CT_proj_from_vc = ps.ProjectAwayFromPowerFixedOr(CT_resh(:), HM.gain, pwr_rnk);
[CS, IND] = ps.PSIICOS_ScanFast(HM.gain, (CT_proj_from_vc));

% CT_rand = rand(size(CT_resh)) *  + j * rand(size(CT_resh));
% [Cs, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_rand), HM.gain, pwr_rnk, Upwr);
% [CS, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_resh), HM.gain, pwr_rnk, Upwr);

con_inds = ups.threshold_connections(CS, threshold, IND);

con = ups.Bundles(con_inds, HM, CtxHHR);
% con_c = con.Clusterize(10,0.02);
figure;
con.Plot(0.2);
