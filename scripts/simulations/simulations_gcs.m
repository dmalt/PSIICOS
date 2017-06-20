% Evaluate performance of PSIICOS compared to GCS %
% ----------------------------------------------- %

% ---- setup simulations ----- %
phi_lag = pi / 2 - pi / 20;
GainSVDTh = 0.01;
n_tr = 100;
ind_scale = 0.5;
evo_scale = 0;
is_usecache = false;
% ---------------------------- %

% ----- setup psiicos ----- %
SL_rnk = 350;
sig_rnk = 0; % corresponds to mean cross-spectrum
cp_part = 'full';
is_fast = true;
Upwr = [];
seed_ind = [];
% ------------------------- %
%
% ------- setup DICS ------ %
lambda = 10;
% ------------------------- %

% ----- setup plotting ---- %
threshold_ps = 300;
threshold_gcs = 300;
% ------------------------- %

[HM, CT, Trials, Ctx] = ups.SimulateData(phi_lag, n_tr, GainSVDTh, ind_scale, evo_scale, is_usecache);

CT_reshape = reshape(mean(CT, 2), sqrt(size(CT,1)), sqrt(size(CT,1)));

[Cs_ps, ~] = ps.PSIICOS_DICS(mean(CT,2), HM.gain, lambda, SL_rnk);

% ------------------- perform computations -------------------- %
[corr, Cp, Upwr]   = ps.PSIICOS(CT, HM.gain, SL_rnk, sig_rnk, Upwr,...
                                  seed_ind, cp_part, is_fast);
[Cs_gcs, ~] = ups.conn.GCS_DICS(CT_reshape, HM.gain, lambda);
[A, Ps, Cs_dics, IND] = ups.conn.DICS(CT_reshape, HM.gain, lambda, true);
% ------------------------------------------------------------- %

IND_dics_gcs = ups.threshold_connections(Cs_gcs, threshold_gcs, IND);
IND_dics_ps  = ups.threshold_connections(Cs_ps, threshold_gcs, IND);
ij_ps  = ups.threshold_connections(corr.data, threshold_ps, corr.IND);



con_ps = ups.Connections(ij_ps, HM, Ctx);
con_dics_gcs = ups.Connections(IND_dics_gcs, HM, Ctx);
con_dics = ups.Connections(IND_dics_ps, HM, Ctx);

% con_ps = con_ps.Clusterize(10, 0.02);
% con_dics_gcs = con_dics_gcs.Clusterize(10, 0.02);
% con_dics = con_dics.Clusterize(10, 0.02);

% plot_dics;
