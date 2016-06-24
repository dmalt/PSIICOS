Rnk = 350;
load GLowRes;
% if exist('Sim.mat', 'file')
% 	sim = load('Sim.mat');
% 	C_long = sim.C_long;
% 	GLR = sim.GLR;
% 	GHR = sim.GHR;
% 	N_ch = sqrt(size(C_long,1));
% 	CT = reshape(sum(C_long, 2), N_ch, N_ch);
% 	GridLocLR = sim.GridLocLR;
% 	GridLocHR = sim.GridLocHR;	
% else
% 	sim = load('InputData4Simulations');
% 	[GLR, C_long, UP] = GenerData(pi / 20, 0.35);
% 	ChUsed = 1:306; ChUsed(3:3:end) = [];
% 	Nsites = size(sim.G.GridLoc, 1);
% 	GainHR_reduced = ReduceToTangentSpace(Nsites, sim.G.Gain, ChUsed);
% 	GHR = UP * GainHR_reduced;
% 	N_ch = sqrt(size(C_long, 1));
% 	CT = reshape(sum(C_long, 2), N_ch, N_ch);
% 	GridLocLR = GLowRes.GridLoc;
% 	GridLocHR = sim.G.GridLoc;
% 	save('Sim.mat', 'GLR', 'GHR', 'C_long', 'GridLocLR', 'GridLocHR', '-v7.3');
% end	
% Cond.CT1 = reshape(sum(Cond.CT, 2), 43, 43);
rel_threshold = 0.6;
[HM, Cp, ~, Ctx] = GenerData(pi / 3, 100);
[INDrap, Cp, Upwr, Cs] = T_PSIICOS(Cp, HM.gain, rel_threshold);
freqBand = [19,23];
timeRange = [0, 0.7];
condName = 'test';
sim_son = Connections('test', INDrap, freqBand, timeRange, Cp, condName, HM, Ctx);

sim_son.PlotCon()
% % INDrap = INDrap(1:10:end,:)
% Ctx = load('/home/dmalt/PSIICOS_osadtchii/anat/0003_pran/tess_cortex_concat.mat');
% drawConnectionsOnBrain(INDrap, R, 1, Ctx);