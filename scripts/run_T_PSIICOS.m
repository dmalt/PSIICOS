rel_threshold = 0.9;
SigRnk = 0;
% size(C_long)
% [INDrap, Cp, Upwr, Cs] = T_PSIICOS(C_long, GLR, rel_threshold, Rnk, SigRnk);
% [indep_topo, c_ss_hat, PVU, SubC, INDrap, Cp, Upwr] = RAP_PSIICOS_Fast(CT, GLR, RAPIts, Rnk);%, Upwr);

% R = GLowRes.GridLoc;
% drawSingleRun(INDrap, R, 'y', Ctx);
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
CurBand = [22, 26];
bInducedOnly = true;
TimeRange = [0.2, 0.8];
real_data = 'Real_data.mat';
% if exist(real_data, 'file')
	% load(real_data)
% else
	Condition = 2;
	trials = LoadTrials(AllSubjects{2}, num2str(Condition), CurBand, TimeRange);
	HM = LoadHeadModel(AllSubjects{2}, num2str(Condition));
	G = HM.gain;
	R = HM.GridLoc;
	CT = CrossSpectralTimeseries(trials.data, bInducedOnly);
	% ERP = mean(trials.data, 3);


% 	[ConData, G] = PrepRealData(AllSubjects{10}, CurBand, bInducedOnly, TimeRange);
% 	R = ConData{Condition}.HM_LR.GridLoc;
	
% 	% save(real_data, 'ConData', 'G', '-v7.3');
% % end
% % [INDrap, Cp, Upwr, Cs] = T_PSIICOS((ConData{Condition}.CrossSpecTime), G, rel_threshold, Rnk, SigRnk);
[INDrap, Cp, Upwr, Cs] = T_PSIICOS((CT), G, rel_threshold, Rnk, SigRnk);

con = Connections(AllSubjects{1}, INDrap, CurBand,...
		                           TimeRange, CT, num2str(Condition));
con.PlotCon()