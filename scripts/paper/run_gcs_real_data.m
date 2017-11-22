AllSubjects = { '0003_pran', ... % 1
                '0019_shev', ... % 2
                '0030_koal', ... % 3
                '0108_bami', ... % 3
                '0062_peek', ... % 5
                '0074_kuni', ... % 6
                '0100_kase', ... % 7
                '0106_supo', ... % 8
                '0109_zvma', ... % 9
                '0130_hagr'};    % 10

subjID = AllSubjects{1};

freqBand = [16, 25];
t_range = [0.4, 0.7];
cond_main = '2';
cond_proj = '1';
GainSVDTh = 0.01;
isInducedOnly = true;
protocolPath = '/home/dmalt/PSIICOS_osadtchii';
isLR = true;
nResamp = 100;
pwr_rnk = 500;
threshold_ps = 100;
threshold_gcs = 100;
lambda = 20000;
SigRnk = 0;
Upwr = [];

HM = ups.LoadHeadModel(subjID, cond_main, protocolPath, isLR, GainSVDTh);

CT_main = ups.GetCTS(subjID, cond_main, freqBand, t_range, GainSVDTh, protocolPath, isInducedOnly, true);
CT_proj = ups.GetCTS(subjID, cond_proj, freqBand, t_range, GainSVDTh, protocolPath, isInducedOnly); 
% CT = ps.ProjFromCond(CT_main, CT_proj);
CT = CT_main - CT_proj;
% CT = CT_main;
CT_reshape = reshape(mean(CT, 2), sqrt(size(CT,1)), sqrt(size(CT,1)));


[Cs_gcs, IND]         = ups.GCS_DICS(CT_reshape, HM.gain, lambda);
[Cs_ps, IND]          = ups.PSIICOS_DICS(CT_reshape(:), HM.gain, lambda, pwr_rnk);
[A, Ps, Cs_dics, IND] = ups.DICS(CT_reshape, HM.gain, lambda);

IND_dics_gcs = ups.threshold_connections(Cs_gcs, threshold_gcs, IND);
IND_dics_ps = ups.threshold_connections(Cs_ps, threshold_ps, IND);
IND_dics = ups.threshold_connections(Cs_dics, threshold_ps, IND);

% [IND_ps,~,Upwr,~] = ps.T_PSIICOS((CT), HM.gain, threshold_ps, pwr_rnk, SigRnk);

[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subjID);

con_ps = ups.Bundles(IND_dics_ps, HM, CtxHHR);
con_dics_gcs = ups.Bundles(IND_dics_gcs, HM, CtxHHR);
con_dics = ups.Bundles(IND_dics, HM, CtxHHR);

plot_dics;
