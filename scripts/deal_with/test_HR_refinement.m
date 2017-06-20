% load ConData1;
% C_long = C.CrossSpecTimeInd;
% % CT = reshape(sum(C.CrossSpecTimeInd(:,90:110), 2), 43, 43);
% % [u, s, v] = svd(C_long);
% % CT = reshape(u(:,1), 43, 43);
% CT = C_long;
% C.G2dHRU = C.UP * C.G2dHR;
% GLR = C.G2dLRU;
% GHR = C.G2dHRU;
% GridLocLR = C.HM_LR.GridLoc;
% GridLocHR = C.HM_HR.GridLoc;
% Upwr = C.Upwr;
path(path, '/home/dmalt/Dropbox/Documents/Education/MEG/Osadchii/gopsiicos/gopsiicos_source/input');
path(path, '/home/dmalt/Dropbox/Documents/Education/MEG/Osadchii/gopsiicos/gopsiicos_source/source/drawscripts');
path(path, '/home/dmalt/Dropbox/Documents/Education/MEG/Osadchii/Simulations/');
RAPIts = 3;
Rnk = 350;
NetworkPairIndex{1} = [1,2];
NetworkPairIndex{2} = [1,2,3];
NPI = NetworkPairIndex{2};
XYZGen = 1.3 * [0.05, 0.04, 0.05; 0.05, -0.04, 0.05; -0.05, 0.04, 0.05; -0.05, -0.04, 0.05; 0.00, 0.05, 0.06; 0.00, -0.05, 0.06];
% iSubj = 1;
load GLowRes;
if exist('Sim.mat', 'file')
	sim = load('Sim.mat');
	C_long = sim.C_long;
	GLR = sim.GLR;
	GHR = sim.GHR;
	% [u, s, v] = svd(C_long);
	N_ch = sqrt(size(C_long,1));
	% CT = reshape(u(:,1), N_ch, N_ch);
	% CT = reshape(sum(C_long, 2), N_ch, N_ch);
	Ctx = sim.Ctx;
	CT = C_long;
	GridLocLR = sim.GridLocLR;
	GridLocHR = sim.GridLocHR;	
	% Gain = sim.GHR.Gain;
	% ChUsed = 1:306; ChUsed(3:3:end) = [];
else
	sim = load('InputData4Simulations');
	[GLR, C_long, UP] = GenerData(pi / 20, 0.35);
	CT = C_long;
	ChUsed = 1:306; ChUsed(3:3:end) = [];
	Nsites = size(sim.G.GridLoc, 1);
	GainHR_reduced = ReduceTangentSpace(Nsites, sim.G.Gain, ChUsed);
	GHR = UP * GainHR_reduced;
	N_ch = sqrt(size(C_long,1));
	CT = reshape(sum(C_long, 2), N_ch, N_ch);
	GridLocLR = GLowRes.GridLoc;
	GridLocHR = sim.G.GridLoc;
	Ctx = sim.Ctx;
	save('Sim.mat', 'GLR', 'GHR', 'C_long', 'GridLocLR', 'GridLocHR', 'Ctx', '-v7.3');
end	
% Cond.CT1 = reshape(sum(Cond.CT, 2), 43, 43);
% [indep_topo, c_ss_hat, PVU, SubC, INDrap, Cp, Upwr] = ...
%  RAP_PSIICOS_Fast(CT, GLR, GHR, GridLocLR, GridLocHR, RAPIts, Rnk); % , Upwr);
% [indep_topo, c_ss_hat, PVU, SubC, INDrap, Cp, Upwr] = RAP_PSIICOS_Fast(CT, GLR, RAPIts, Rnk);%, Upwr);
% [indep_topo, c_ss_hat, PVU, SubC, INDrap, Cp, Upwr] = T_PSIICOS(CT, GLR, Rnk);
Dmax = 0.015;
N = 100;
Nsites = size(GLR, 2) / 2;
R = GLowRes.GridLoc;
IND = UpperDiagToPairs(Nsites);
Labels = {'k.-','b.-','ro-','bo-','m.'};
CT = C_long;
[INDrap, Cp, Upwr, Cs] = T_PSIICOS_subcorr(CT, GLR, Rnk, 3, Upwr);
Q = Cs;
[SPC, TPR] = GenerateROC(Q, Dmax, R, IND, N, XYZGen, NPI); 
plot(1 - SPC, TPR, 'LineWidth', 2);
hold on;

for SignalRnk = 1:5
	[INDrap, Cp, Upwr, Cs] = T_PSIICOS(CT, GLR, Rnk, SignalRnk, Upwr);
	Q = Cs;
	[SPC, TPR] = GenerateROC(Q, Dmax, R, IND, N, XYZGen, NPI); 
	plot(1 - SPC, TPR, 'LineWidth', 2);
	hold on;
end
CT = sum(C_long, 2);
[INDrap, Cp, Upwr, Cs] = T_PSIICOS(CT, GLR, Rnk, 1, Upwr);
Q = Cs;
[SPC, TPR] = GenerateROC(Q, Dmax, R, IND, N, XYZGen, NPI); 
plot(1 - SPC, TPR, 'LineWidth', 2);
% hold on;

legend('subcorr',num2str(1), num2str(2), num2str(3), num2str(4), num2str(5), 'sum');
% % legend(num2str(1), num2str(2), num2str(3), num2str(4), num2str(5), num2str(6), num2str(7), num2str(8), num2str(9), num2str(10), 'sum');
% % legend(num2str(1), num2str(2), num2str(3), num2str(4), num2str(5));%, num2str(6), num2str(7), num2str(8), num2str(9), num2str(10));
% % legend(num2str(1), num2str(3), num2str(5), num2str(7), num2str(9));
% axis([0 0.003 0 1.1])

		
% drawSingleRun(INDrap, R, 'y', Ctx);