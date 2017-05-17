
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
% ---- setup psiicos ----- %
SL_rnk = 350;
sig_rnk = 0;
cp_part = 'real';
is_fast = false;
Upwr = [];
% ------------------------ %

% ------------- load individual and avg brain surfaces -------------------- %
[Ctx_dst, CtxHR_dst, CtxHHR_dst] = ups.GetCtx('@default_subject', protocol_path);
n_src_dst = length(Ctx_dst.Vertices);

[Ctx_src, CtxHR_src, CtxHHR_src] = ups.GetCtx(subj_ID, protocol_path);
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
% ---------------------------------------------- %

threshold = 100;


HM = ups.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);

[~, ~, Ps_pre_dst, CS_pre_dst] = prepare_cond(subj_ID, condition, protocol_path,...
                                              HM, GainSVDTh, freq_band,...
                                              time_range_pre, lambda, seed_src_ind,...
                                              SL_rnk, sig_rnk, Upwr, cp_part,...
                                              is_fast, Wmat);

[~, ~, Ps_post_dst, CS_post_dst] = prepare_cond(subj_ID, condition, protocol_path,...
                                                HM, GainSVDTh, freq_band,...
                                                time_range_post, lambda,...
                                                seed_src_ind, SL_rnk, sig_rnk,...
                                                Upwr, cp_part, is_fast, Wmat);


CS_dst = (CS_post_dst.data - CS_pre_dst.data) ./ CS_pre_dst.data;
% CS_dst = (CS_post_dst - CS_pre_dst) ./ CS_pre_dst;


% figure;
% h = plot_brain_cmap(CtxHHR, Ctx, [], Ps_pre, zeros(size(Ps_pre)), 0.2);

figure;
h = plot_brain_cmap_hemisplit(CtxHHR_dst, Ctx_dst, [], CS_dst,...
                              zeros(size(CS_dst)), 0, seed_dst_xyz_approx);
% CT_proj_from_vc = ps.ProjectAwayFromPowerFixedOr(CT_resh(:), HM.gain, pwr_rnk);
% [CS, IND] = ps.PSIICOS_ScanFast(HM.gain, (CT_proj_from_vc));

% CT_rand = rand(size(CT_resh)) *  + j * rand(size(CT_resh));
% [Cs, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_rand), HM.gain, pwr_rnk, Upwr);
% [CS, IND, Cp, Upwr] = ps.PSIICOS_Fast_Norm((CT_resh), HM.gain, pwr_rnk, Upwr);

% con_inds = ups.threshold_connections(CS, threshold, IND);

% con = ups.Connections(subj_ID, con_inds, freq_band, time_range, CT, condition, HM, CtxHHR);
% con_c = con.Clusterize(10,0.02);
% figure;
% con.Plot(0.2);
