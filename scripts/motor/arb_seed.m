protocol_path = '/home/dmalt/fif_matlab/Brainstorm_db/motorprobe';
subj_ID = 'AB';
% Found sma for ab, rv, nv
suffix = '_control_RH_raw_tsss_ica';

condition = [subj_ID, suffix];
freq_band = [19, 25];
time_range_pre = [-0.6, -0.2];
time_range_post = [0, 0.4];

% n_resamp = 2000;
n_resamp = 1;
cond_proj_rnk = 100;

% ---- setup head model loading params ---- %
GainSVDTh = 0.0001;
isLR = true;
% ----------------------------------------- %

lambda = 1;

% ---- setup psiicos ----- %
SL_rnk = 1000;
sig_rnk = 20;
cp_part = 'real';
is_fast = false;
Upwr = [];
% ------------------------ %

% --------------- load individual and avg brain surfaces -------------------- %
[Ctx_dst, CtxHR_dst, CtxHHR_dst] = ups.bst.GetCtx('@default_subject', protocol_path);
n_src_dst = length(Ctx_dst.Vertices);

[Ctx_src, CtxHR_src, CtxHHR_src] = ups.bst.GetCtx(subj_ID, protocol_path);
n_src_src = length(Ctx_src.Vertices);

Wmat = map_on_default(Ctx_src, Ctx_dst);
% ---------------------------------------------------------- %

% ------ set seed location on avg brain -------- %
% Uncomment one of these lines to test for high SNR effect
l_M1_dst_xyz = [-0.002, 0.0351, 0.1003]; % good seed
% l_M1_dst_xyz = [0.102, 0.0351, 0.1003]; % bad seed
% l_M1_dst_xyz = [-0.102, 0.0351, 0.1003]; % bad seed

seed_dst_xyz_approx = l_M1_dst_xyz;
seed_dst_ind = ups.FindXYZonGrid(seed_dst_xyz_approx, Ctx_dst.Vertices);
seed_dst_xyz = Ctx_dst.Vertices(seed_dst_ind,:);
seed_dst_indicator = zeros(n_src_dst, 1);
seed_dst_indicator(seed_dst_ind) = 1;
[~, seed_src_ind]= max(Wmat' * seed_dst_indicator);
seed_src_xyz = Ctx_src.Vertices(seed_src_ind,:);
% ---------------------------------------------- %

threshold = 100;


HM = ups.bst.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);


[tr_pre, tr_post,...
 CT_post, Ps_post_dst,...
 CS_post_dst_proj, Upwr, ~] = prepare_cond(subj_ID, condition,...
                                        protocol_path,...
                                        HM, freq_band,...
                                        time_range_pre, time_range_post,...
                                        lambda,...
                                        seed_src_ind, SL_rnk, sig_rnk,...
                                        Upwr, cp_part, is_fast, Wmat, seed_dst_ind);

% CS_dst = (CS_post_dst.data - CS_pre_dst.data);
CS_dst = CS_post_dst_proj;


% clust_threshold = 0.12;
clust_threshold = 0.7;

CS_dst_perm = zeros(n_src_dst, n_resamp);
adj_mat = get_adjacency_matrix(Ctx_dst);
clust_stat = cell(n_resamp,1);





% CS_dst = (CS_post_dst - CS_pre_dst) ./ CS_pre_dst;

thresh_mask = CS_dst.data < 0.9;
figure;
h = plot_brain_cmap_hemisplit(CtxHHR_dst, Ctx_dst, [], CS_dst.data,...
                              zeros(size(CS_dst.data)), 0.1, seed_dst_xyz_approx);
