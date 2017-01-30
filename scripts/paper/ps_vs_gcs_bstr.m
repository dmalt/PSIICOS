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
nResamp = 300;
pwr_rnk = 500;
threshold = 50;
nTrials_used = 10;
Upwr = []; 
bInducedOnly = true;
threshold_gcs = 50; 

HM = ups.LoadHeadModel(subjID, cond_main, protocolPath, isLR, GainSVDTh);
trials = ups.LoadTrials(subjID, cond_main, freqBand, t_range, GainSVDTh, protocolPath);
CtxHR = load('/home/dmalt/PSIICOS_osadtchii/anat/@default_subject/tess_cortex_pial_high.mat');
save_fname = ['boots_CT_cache', num2str(nResamp), '_', num2str(bInducedOnly), '_', num2str(nTrials_used), '_', cond_main, '.mat'];

if ~exist(save_fname)
	[boots_CT, resamples] = ups.bootstrap_CT(trials.data, nResamp, bInducedOnly, nTrials_used);
	save(save_fname, 'boots_CT', 'resamples', '-v7.3');
else
	load(save_fname);
end

lambdas = [0.001, 0.01, 0.1, 1, 10, 40, 100, 1000];
iLambda = 1;

for lambda = lambdas
	% [Ctx, CtxHR]= ups.GetCtx(subjID);
	% boots_IND = ups.bootstrap_PSIICOS(trials.data, HM.gain, nResamp, pwr_rnk, threshold, 'full', isInducedOnly);


	boots_IND_gcs = cell(1, nResamp);
	boots_IND_ps  = cell(1, nResamp);

	for iCT = 1:length(boots_CT)
		iCT
		tic
		fprintf('ps')
		[boots_IND_ps{iCT}, ~, Upwr, ~] = ps.T_PSIICOS(mean(boots_CT{iCT}, 2), HM.gain, threshold, pwr_rnk, 0, Upwr);
		toc
		CT_reshape = reshape(mean(boots_CT{iCT}, 2), sqrt(size(boots_CT{iCT},1)), sqrt(size(boots_CT{iCT},1)));
		tic
		fprintf('gcs')
		[Cs_gcs, IND]             = ups.GCS_DICS(CT_reshape, HM.gain, lambda);
		toc
		boots_IND_gcs{iCT} = ups.threshold_connections(Cs_gcs, threshold_gcs, IND);
	end

	con_ps_all = ups.Connections(subjID, boots_IND_ps, freqBand, t_range, [], '2', HM, CtxHR);
	con_clust_av_ps = con_ps_all.ClustAndAvCells(1,0.02);
	con_m_ps = con_clust_av_ps.Merge();
	con_ps{iLambda} = con_m_ps.Clusterize(20, 0.013);

	% con_c_ps.PlotViews(0.2, 2, 0.003);

	con_gcs_all = ups.Connections(subjID, boots_IND_gcs, freqBand, t_range, [], '2', HM, CtxHR);
	con_clust_av_gcs = con_gcs_all.ClustAndAvCells(1,0.02);
	con_m_gcs = con_clust_av_gcs.Merge();
	con_dics_gcs{iLambda} = con_m_gcs.Clusterize(20, 0.013);
	iLambda = iLambda + 1;
end


plot_gcs
% con_c_gcs.PlotViews(0.2, 2, 0.003);
% con = ups.Connections(subjID, boots_IND, freqBand, t_range, [], '2', HM, CtxHR);
% con_clust_av = con.ClustAndAvCells(1,0.02);
% con_m = con_clust_av.Merge();
% con_c = con_m.Clusterize(20, 0.011);
% con_c.PlotViews(0.2, 2, 0.003);
% CT_main = ups.GetCTS(subjID, cond_main, freqBand, t_range, GainSVDTh, protocolPath, isInducedOnly);
% CT_proj = ups.GetCTS(subjID, cond_proj, freqBand, t_range, GainSVDTh, protocolPath, isInducedOnly);
% CT_final = ps.ProjFromCond(CT_main, CT_proj, 6);
% CT_final = CT_main;
% [INDrap, Cp, Upwr, corr] = ps.T_PSIICOS(CT_final, HM.gain, threshold, pwr_rnk, 0, Upwr);
% con = ups.Connections(subjID, INDrap, freqBand, [], [], '2', HM, CtxHR);
