clear all;
subjID = 'test';
PhaseLag = pi / 20;

[HM, CT, Trials, Ctx] = ups.SimulateData(PhaseLag, 100, 0.01, 0.5, 0, false);

freqBand = [8, 12];
t_range = [0.4, 0.7];
GainSVDTh = 0.01;
isInducedOnly = true;
isLR = true;
pwr_rnk = 500;
threshold_ps = 300;
threshold_gcs = 300;
SigRnk = 0;
lambda = 100;
Upwr = [];

CT_reshape = reshape(mean(CT, 2), sqrt(size(CT,1)), sqrt(size(CT,1)));

[Cs_ps, IND]   = ups.PSIICOS_DICS(mean(CT,2), HM.gain, lambda, pwr_rnk);
[Cs_gcs, IND] = ups.GCS_DICS(CT_reshape, HM.gain, lambda);
[A, Ps, Cs_dics, IND] = ups.DICS(CT_reshape, HM.gain, lambda);

IND_dics_gcs = ups.threshold_connections(Cs_gcs, threshold_gcs, IND);
IND_dics_ps  = ups.threshold_connections(Cs_ps, threshold_gcs, IND);
IND_dics  = ups.threshold_connections(Cs_dics, threshold_gcs, IND);

% [IND_ps,~,Upwr,~] = ps.T_PSIICOS((CT), HM.gain, threshold_ps, pwr_rnk, SigRnk);


con_ps       = ups.Connections(subjID, IND_dics_ps,  freqBand, t_range, CT, 'test', HM, Ctx);
con_dics_gcs = ups.Connections(subjID, IND_dics_gcs, freqBand, t_range, CT, 'test', HM, Ctx);
con_dics = ups.Connections(subjID, IND_dics, freqBand, t_range, CT, 'test', HM, Ctx);

plot_dics;
% con_dics_gcs.Plot();
