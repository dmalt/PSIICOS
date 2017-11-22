AllSubjects = { '0003_pran', ... % 1
				'0019_shev', ... % 2
				'0030_koal', ... % 3
				'0108_bami', ... % 4
				'0062_peek', ... % 5
				'0074_kuni', ... % 6
				'0100_kase', ... % 7
				'0106_supo', ... % 8
				'0109_zvma', ... % 9
				'0130_hagr'};    % 10

subjID = AllSubjects{5};

freqBand = [16, 25];
t_range = [0., 0.6];
cond_main = '4';
cond_proj = '1';
GainSVDTh = 0.01;
isInducedOnly = true;
protocolPath = '/home/dmalt/PSIICOS_osadtchii';
isLR = true;
nResamp = 100;
pwr_rnk = 500;
threshold = 1000;
SigRnk = 0;
Upwr = [];

HM = ups.LoadHeadModel(subjID, cond_main, protocolPath, isLR, GainSVDTh);

CT_main = ups.GetCTS(subjID, cond_main, freqBand,...
                     t_range, GainSVDTh, protocolPath,...
                     isInducedOnly, true);

CT_proj = ups.GetCTS(subjID, cond_proj, freqBand,...
                     t_range, GainSVDTh, protocolPath,...
                     isInducedOnly); 

% CT = ps.ProjFromCond(CT_main, CT_proj);
norm(CT_main)
norm(CT_proj)
CT = CT_main - CT_proj;
norm(CT)
CT_reshape = reshape(mean(CT, 2), sqrt(size(CT, 1)), sqrt(size(CT, 1)));

% CtxHR = load('/home/dmalt/PSIICOS_osadtchii/anat/@default_subject/tess_cortex_pial_high.mat');
[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subjID);
% [IND_music_tot,~,Upwr,~] = ps.T_PSIICOS((CT), HM.gain, threshold, pwr_rnk, SigRnk);
[CS, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_reshape), HM.gain, pwr_rnk, Upwr);

con_inds = ups.threshold_connections(CS, threshold, IND);

con = ups.Bundles(con_inds, HM, CtxHHR);
% con_c = con.Clusterize(10,0.02);
figure;
con.Plot(0.2);
