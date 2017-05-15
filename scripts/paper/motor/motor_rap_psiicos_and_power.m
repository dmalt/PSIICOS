
protocol_path = '/home/dmalt/fif_matlab/Brainstorm_db/motorprobe';

AllSubjects = { 'AP', ... % 1
				'MG', ... % 2
				'MT', ... % 3
				'NV', ... % 4
				'RV', ... % 5
				'SR', ... % 6
				'VK', ... % 7
				'VM'};    % 8
% subj_ID = 'AP';
suffix = '_control_RH_raw_tsss_ica';

% condition = [subj_ID, suffix]; 
freq_band = [10, 20];
time_range = [-1, -0.6];
GainSVDTh = 0.01;
isLR = true;
lambda = 30;
pwr_rnk = 350;
isInducedOnly = true;
Upwr = [];
threshold = 100;
RAPIts = 3;

% for i_subj = 1:length(AllSubjects)

    subj_ID = AllSubjects{6};
    condition = [subj_ID, suffix]; 

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


    [indep_topo,...
     c_ss_hat,...
     PVU,...
     SubC,...
     INDrap_true,...
     Cp,...
     Upwr,...
     Cs_true,...
     return_qp] = ps.RAP_PSIICOS_Fast(mean(CT, 2), HM.gain, RAPIts, pwr_rnk);

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

% end
