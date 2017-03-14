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

% These parameters give nice r_motor - l_motor, r_motor - r_audio network
% subjID = AllSubjects{9};
% freqBand = [8, 12];
% t_range = [0.4, 0.7];

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
threshold = 100;
SigRnk = 0;
Upwr = [];

HM = ups.LoadHeadModel(subjID, cond_main, protocolPath, isLR, GainSVDTh);

CT_main = ups.GetCTS(subjID, cond_main, freqBand, t_range,...
                     GainSVDTh, protocolPath, isInducedOnly, false);
% CT_proj = ups.GetCTS(subjID, cond_proj, freqBand, t_range, GainSVDTh, protocolPath, isInducedOnly); 
% CT = ps.ProjFromCond(CT_main, CT_proj);
CT = CT_main;

% CtxHR = load('/home/dmalt/PSIICOS_osadtchii/anat/@default_subject/tess_cortex_pial_high.mat');
[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subjID);

[IND_music_tot, Cpvec, Upwr, corr, IND_run] = ps.T_PSIICOS(CT, HM.gain,...
                                                           threshold, pwr_rnk, SigRnk, Upwr);
con = ups.Connections(subjID, IND_music_tot, freqBand,...
                      t_range, CT, cond_main, HM, CtxHHR);
% con_c = con.Clusterize(10,0.02);
figure;
con.Plot(0.2);
