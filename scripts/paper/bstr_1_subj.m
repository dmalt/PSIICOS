subjID = '0003_pran';
% subjID = '0100_kase';
% subjID = '0019_shev';
% subjID = '0109_zvma';


freqBand = [16,25];
t_range = [0.4, 0.7];
cond_main = '2';
cond_proj = '1';
GainSVDTh = 0.01;
isInducedOnly = true;
protocolPath = '/home/dmalt/PSIICOS_osadtchii';
isLR = true;
nResamp = 500;
pwr_rnk = 500;
threshold = 0.97;
Upwr = []; 

trials = ups.LoadTrials(subjID, cond_main, freqBand, t_range, GainSVDTh, protocolPath);
HM = ups.LoadHeadModel(subjID, cond_main, protocolPath, isLR, GainSVDTh);
% [Ctx, CtxHR]= ups.GetCtx(subjID);
CtxHR = load('/home/dmalt/PSIICOS_osadtchii/anat/@default_subject/tess_cortex_pial_high.mat');
boots_IND = ups.bootstrap_PSIICOS(trials.data, HM.gain, nResamp, pwr_rnk, threshold, 'full', isInducedOnly);
con = ups.Bundles(boots_IND, HM, CtxHR);
con_clust_av = con.ClustAndAvCells(1,0.02);
con_m = con_clust_av.Merge();
con_c = con_m.Clusterize(20, 0.011);
con_c.PlotViews(0.2, 2, 0.003);
% CT_main = ups.GetCTS(subjID, cond_main, freqBand, t_range, GainSVDTh, protocolPath, isInducedOnly);
% CT_proj = ups.GetCTS(subjID, cond_proj, freqBand, t_range, GainSVDTh, protocolPath, isInducedOnly);
% CT_final = ps.ProjFromCond(CT_main, CT_proj, 6);
% CT_final = CT_main;
% [INDrap, Cp, Upwr, corr] = ps.T_PSIICOS(CT_final, HM.gain, threshold, pwr_rnk, 0, Upwr);
% con = ups.Bundles(INDrap, HM, CtxHR);
