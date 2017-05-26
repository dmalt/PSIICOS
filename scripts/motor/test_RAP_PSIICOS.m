protocol_path = '/home/dmalt/fif_matlab/Brainstorm_db/motorprobe';
subj_ID = 'AB';
suffix = '_control_RH_raw_tsss_ica';

condition = [subj_ID, suffix];
freq_band = [19, 25];
time_range_pre = [-0.6, -0.2];
time_range_post = [0, 0.4];

GainSVDTh = 0.01;
isLR = true;
lambda = 1;
% ---- setup rap_psiicos ----- %
SL_rnk = 350;
sig_rnk = 20;
cp_part = 'full';
is_fast = false;
Upwr = [];
n_rap = 5;
% ------------------------ %


% ------------- load individual and avg brain surfaces -------------------- %
[Ctx_dst, CtxHR_dst, CtxHHR_dst] = ups.bst.GetCtx('@default_subject', protocol_path);
n_src_dst = length(Ctx_dst.Vertices);

[Ctx_src, CtxHR_src, CtxHHR_src] = ups.bst.GetCtx(subj_ID, protocol_path);
n_src_src = length(Ctx_src.Vertices);

Wmat = map_on_default(Ctx_src, Ctx_dst);
% ---------------------------------------------------------- %

% ------ set seed location on avg brain -------- %
l_M1_dst_xyz = [-0.002, 0.0351, 0.1003];
seed_dst_xyz_approx = l_M1_dst_xyz;
seed_dst_ind = ups.FindXYZonGrid(seed_dst_xyz_approx, Ctx_dst.Vertices);
seed_dst_xyz = Ctx_dst.Vertices(seed_dst_ind,:);
seed_dst_indicator = zeros(n_src_dst, 1);
seed_dst_indicator(seed_dst_ind) = 1;
[~, seed_src_ind]= max(Wmat' * seed_dst_indicator);
seed_src_xyz = Ctx_src.Vertices(seed_src_ind,:);


HM = ups.bst.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);

trials_pre = ups.bst.LoadTrials(subj_ID, condition,...
                            freq_band, time_range_pre,...
                            HM, protocol_path);

trials_post = ups.bst.LoadTrials(subj_ID, condition,...
                            freq_band, time_range_post,...
                            HM, protocol_path);

tr_pre = trials_pre.data;
tr_post = trials_post.data;
CT_pre = ups.conn.CrossSpectralTimeseries(tr_pre, true);
CT_post = ups.conn.CrossSpectralTimeseries(tr_post, true);


% corr_pre = ps.RAP_PSIICOS(CT_pre, HM.gain, SL_rnk,...
%                       sig_rnk, Upwr, seed_src_ind,...
%                       cp_part, is_fast, n_rap);

CT_post = ps.ProjFromCond(CT_post, CT_pre, 9);

corr_post = ps.RAP_PSIICOS(CT_post, HM.gain, SL_rnk,...
                      sig_rnk, Upwr, seed_src_ind,...
                      cp_part, is_fast, n_rap);

% for i = 1:n_rap
%     corr{i} = Wmat * ((corr_post{i}.data - corr_pre{i}.data));
% end

% check_bad = Ctx_dst.Vertices(1315,:);

figure;
h = plot_brain_cmap_hemisplit(CtxHHR_dst, Ctx_dst, [], Wmat * corr_post{1}.data,...
                              zeros(size(corr{1})), 0.1, seed_dst_xyz);
