clear all;
subjID = 'test';
% PhaseLag = pi / 20;
PhaseLag = pi / 2 - pi / 20;

[HM, CT, Trials, Ctx] = ups.SimulateData(PhaseLag, 100, 0.01, 1, 0, false);

freqBand = [2, 20];
t_range = [0.4, 0.7];
GainSVDTh = 0.01;
isInducedOnly = true;
isLR = true;
pwr_rnk = 350;
threshold_ps = 300;
threshold_gcs = 300;
SigRnk = 0;
lambda = 1000;
Upwr = [];

CT_reshape = reshape(mean(CT, 2), sqrt(size(CT,1)), sqrt(size(CT,1)));

% [Cs_ps, IND]   = ups.PSIICOS_DICS(mean(CT,2), HM.gain, lambda, pwr_rnk);
[IND_dics_ps, Cs_ps, IND]   = ps.T_PSIICOS(mean(CT,2), HM.gain, threshold_ps, pwr_rnk, 0);
[Cs_gcs, IND] = ups.GCS_DICS(CT_reshape, HM.gain, lambda);
[A, Ps, Cs_dics, IND] = ups.DICS(CT_reshape, HM.gain, lambda, true);

IND_dics_gcs = ups.threshold_connections(Cs_gcs, threshold_gcs, IND);
% IND_dics_ps  = ups.threshold_connections(Cs_ps, threshold_gcs, IND);
IND_dics  = ups.threshold_connections(Cs_dics, threshold_gcs, IND);

% [IND_ps,~,Upwr,~] = ps.T_PSIICOS((CT), HM.gain, threshold_ps, pwr_rnk, SigRnk);


con_ps       = ups.Connections(subjID, IND_dics_ps,  freqBand, t_range, CT, 'test', HM, Ctx);
con_dics_gcs = ups.Connections(subjID, IND_dics_gcs, freqBand, t_range, CT, 'test', HM, Ctx);
con_dics = ups.Connections(subjID, IND_dics, freqBand, t_range, CT, 'test', HM, Ctx);

con_ps = con_ps.Clusterize(10, 0.02);
con_dics_gcs = con_dics_gcs.Clusterize(10, 0.02);
con_dics = con_dics.Clusterize(10, 0.02);

plot_dics;
% con_dics_gcs.Plot();
