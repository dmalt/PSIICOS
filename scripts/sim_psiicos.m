% Test PSIICOS performance on simulated data %
% ------------------------------------------ %

subjID = 'test';
% PhaseLag = pi / 20;
PhaseLag = pi / 2 - pi / 20;
GainSVDTh = 0.01;

[HM, CT, Trials, Ctx] = ups.SimulateData(PhaseLag, 100, GainSVDTh, 1, 0, false);

freqBand = [2, 20];
t_range = [0.4, 0.7];
isInducedOnly = true;
isLR = true;
% --- setup psiicos params --- %
SL_rnk = 350;
sig_rnk = 0; % corresponds to mean cross-spectrum
cp_part = 'full';
is_fast = false;
seed_ind = []; % if empty compute all-to-all
Upwr = []; % if empty recompute projection from SL
% ---------------------------- %
threshold_ps = 300;
threshold_gcs = 300;

CT_reshape = reshape(mean(CT, 2), sqrt(size(CT,1)), sqrt(size(CT,1)));

% -------- compute network subspace correlations with psiicos -------- %
% profile on
tic;
[corr, Cp, Upwr]   = ps.PSIICOS(CT, HM.gain, SL_rnk, sig_rnk, Upwr,...
                                seed_ind, cp_part, is_fast);
toc
% profile viewer;
% -------------------------------------------------------------------- %


ij_ps  = ups.threshold_connections(corr.data, threshold_ps, corr.IND);

con_ps = ups.Connections(ij_ps, HM, Ctx);
con_ps_clust = con_ps.Clusterize(10, 0.02);
con_ps_av = con_ps_clust.Average();
con_ps_av.Plot(0.2, 4, 0.004);