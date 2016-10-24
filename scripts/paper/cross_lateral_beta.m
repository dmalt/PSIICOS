rel_threshold = 0.8;
SigRnk = 0;
Rnk = 350;

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

freqBand = [16, 25];
isInducedOnly = true;
timeRange = [0, 0.7];

subjID = AllSubjects{4};
GainSVDTh = 0.01;
real_data = 'Real_data.mat';
protocolPath = '/home/dmalt/PSIICOS_osadtchii';
condition = 2;
	% trials = ups.LoadTrials(AllSubjects{2}, num2str(condition), freqBand, timeRange);
HM = ups.LoadHeadModel(subjID, num2str(condition));
G = HM.gain;
R = HM.GridLoc;

CT = ups.GetCTS(subjID, num2str(condition), freqBand, timeRange, GainSVDTh, protocolPath, isInducedOnly);
[~,CtxHR] = ups.GetCtx(subjID);

[IND_music_tot,~,Upwr,~] = ps.T_PSIICOS((CT), G, rel_threshold, Rnk, SigRnk);
IND_music_real = ps.T_PSIICOS(real(CT), G, rel_threshold, Rnk, SigRnk, Upwr);
IND_music_imag = ps.T_PSIICOS(imag(CT), G, rel_threshold, Rnk, SigRnk, Upwr);


% CtxHR = load('/home/dmalt/PSIICOS_osadtchii/anat/0019_shev/tess_cortex_pial_low_fig.mat')

con = ups.Connections(subjID, IND_music_tot, freqBand,...
		                           timeRange, CT, num2str(condition), HM, CtxHR);

% --------------------- Plotting part --------------------------- %
subplot(1,3,1)
con_clust = con.Clusterize(20, 0.02);
con_clust.Plot(0.5, 1, 0.004)

con_real = ups.Connections(subjID, IND_music_real, freqBand,...
		                           timeRange, CT, num2str(condition), HM, CtxHR);

subplot(1,3,2)
con_clust_real = con_real.Clusterize(20, 0.02);
con_clust_real.Plot(0.5, 1, 0.004)

con_imag = ups.Connections(subjID, IND_music_imag, freqBand,...
		                           timeRange, CT, num2str(condition), HM, CtxHR);
subplot(1,3,3)
con_clust_imag = con_imag.Clusterize(20, 0.02);
con_clust_imag.Plot(0.5, 1, 0.004)
% ----------------------------------------------------------------- %