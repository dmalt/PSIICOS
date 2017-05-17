% clear all;
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

for i_subj  = 5:10

subj_ID = AllSubjects{i_subj};

protocol_path = '/home/dmalt/PSIICOS_osadtchii';

condition = '2';

freq_band = [16, 25];
time_range = [0.4, 0.7];

GainSVDTh = 0.01;
isLR = true;
lambda = 30;
pwr_rnk = 500;
Upwr = [];
threshold = 100;
isInducedOnly = true;
RAPIts = 10;

trials = ups.LoadTrials(subj_ID, condition,...
                        freq_band, time_range,...
                        GainSVDTh, protocol_path);
trials_1 = ups.LoadTrials(subj_ID, '1',...
                        freq_band, time_range,...
                        GainSVDTh, protocol_path);

HM = ups.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);

tr = trials.data;%(:, :, 51:100);% + 50 * rand(size(trials.data));
tr1 = trials_1.data;
CT = ups.CrossSpectralTimeseries(tr, isInducedOnly);
CT1 = ups.CrossSpectralTimeseries(tr1, isInducedOnly);

CT = CT - CT1;

CT_resh = reshape(mean(CT,2), sqrt(size(CT,1)), sqrt(size(CT,1)));
[A, Ps] = ups.DICS(CT_resh, HM.gain, lambda);
[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subj_ID, protocol_path);

figure;
h = plot_brain_cmap(CtxHHR, Ctx, [], Ps, zeros(size(Ps)), 0.2);


[indep_topo,...
 c_ss_hat,...
 PVU,...
 SubC,...
 INDrap_true,...
 Cp,...
 Upwr,...
 Cs_true,...
 return_qp] = ps.RAP_PSIICOS_Fast(mean(CT, 2), HM.gain, RAPIts, pwr_rnk, Upwr);

con_ps = ups.Connections(subj_ID,...
                              INDrap_true,...
                              freq_band,...
                              time_range,...
                              CT,...
                              'test',...
                              HM,...
                              CtxHHR);
% CT_rand = rand(size(CT_resh)) *  + j * rand(size(CT_resh));
% [Cs, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_rand), HM.gain, pwr_rnk, Upwr);
% [CS, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_resh), HM.gain, pwr_rnk, Upwr);

% con_c = con.Clusterize(10,0.02);
figure;
con_ps.Plot(0.2);
end
