% subjID = '0003_pran';
% subjID = '0100_kase';
% subjID = '0019_shev';
subjNames = {   '0003_pran', ... % 1
				'0019_shev', ... % 2
				'0030_koal', ... % 3
				'0062_peek', ... % 4
				'0074_kuni', ... % 5
				'0100_kase', ... % 6
				'0106_supo', ... % 7
				'0108_bami', ... % 8
				'0109_zvma', ... % 9
				'0130_hagr'};    % 10
subjID = subjNames{1};


freqBand = [16,25];
t_range = [-0.5, 0.];
cond_main = '2';
cond_proj = '1';
GainSVDTh = 0.01;
isInducedOnly = true;
protocolPath = '/home/dmalt/PSIICOS_osadtchii';
isLR = true;
nResamp = 300;
pwr_rnk = 500;
threshold = 50;
nTrials_used = 10;
Upwr = [];
bInducedOnly = true;
threshold_gcs = 50;
lambda = 30;

HM = ups.LoadHeadModel(subjID, cond_main, protocolPath, isLR, GainSVDTh);
trials = ups.LoadTrials(subjID, cond_main, freqBand, t_range, GainSVDTh, protocolPath);
% CtxHR = load('/home/dmalt/PSIICOS_osadtchii/anat/@default_subject/tess_cortex_pial_high.mat');

CT = ups.CrossSpectralTimeseries(trials.data, true);
CT_resh = reshape(mean(CT,2), sqrt(size(CT,1)), sqrt(size(CT,1)));
[A, Ps] = ups.DICS(CT_resh, HM.gain, lambda);
[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subjID);

figure;
h = plot_brain_cmap(CtxHHR, Ctx, [], Ps, zeros(size(Ps)), 0.2);
