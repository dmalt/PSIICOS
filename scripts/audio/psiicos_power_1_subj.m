
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

subj_ID = AllSubjects{9};

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

trials = ups.LoadTrials(subj_ID, condition,...
                        freq_band, time_range,...
                        GainSVDTh, protocol_path);

HM = ups.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);

tr = trials.data;%(:, :, 51:100);% + 50 * rand(size(trials.data));
CT = ups.CrossSpectralTimeseries(tr, isInducedOnly);

CT_resh = reshape(mean(CT,2), sqrt(size(CT,1)), sqrt(size(CT,1)));
[A, Ps] = ups.DICS(CT_resh, HM.gain, lambda);
[Ctx, CtxHR, CtxHHR] = ups.GetCtx(subj_ID, protocol_path);

figure;
h = plot_brain_cmap(CtxHHR, Ctx, [], Ps, zeros(size(Ps)), 0.2);

CT_proj_from_vc = ps.ProjectFromSlComplete(CT_resh(:), HM.gain, pwr_rnk);
[CS, IND] = ps.PSIICOS_ScanFast(HM.gain, (CT_proj_from_vc));

% CT_rand = rand(size(CT_resh)) *  + j * rand(size(CT_resh));
% [Cs, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_rand), HM.gain, pwr_rnk, Upwr);
% [CS, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_resh), HM.gain, pwr_rnk, Upwr);

con_inds = ups.threshold_connections(CS, threshold, IND);

con = ups.Connections(subj_ID, con_inds, freq_band, time_range, CT, condition, HM, CtxHHR);
con_c = con.Clusterize(10,0.02);
figure;
con_c.Plot(0.2);
