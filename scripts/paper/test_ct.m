clear all;
% subjID = '0003_pran';
% subjID = '0100_kase';
% subjID = '0019_shev';
subjID = '0109_zvma';


freqBand = [16,25];
t_range = [0.4, 0.7];
cond_main = '2';
cond_proj = '1';
GainSVDTh = 0.01;
isInducedOnly = true;
protocol_path = '/home/dmalt/PSIICOS_osadtchii';
isLR = true;
nResamp = 300;
pwr_rnk = 350;
threshold = 50;
nTrials_used = 10;
Upwr = [];
bInducedOnly = true;
threshold_gcs = 50;

HM = ups.LoadHeadModel(subjID, cond_main, protocol_path, isLR, GainSVDTh);
trials = ups.LoadTrials(subjID, cond_main, freqBand, t_range, GainSVDTh, protocol_path);
% trials = ups.LoadTrials(subjID, cond_main, freqBand, t_range, GainSVDTh, protocol_path);
% CtxHR = load('/home/dmalt/PSIICOS_osadtchii/anat/@default_subject/tess_cortex_pial_high.mat');
[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subjID, protocol_path);
[~, ~, CtxHHR] = ups.GetCtx(subjID);
tr = trials.data(:,:,110:120);

ct = ups.CrossSpectralTimeseries(tr, false);
con_inds = ps.T_PSIICOS(ct, HM.gain, 100);
con = ups.Bundles(con_inds, HM, CtxHR);
con.Plot();
