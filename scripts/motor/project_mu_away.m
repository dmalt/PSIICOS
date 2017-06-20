
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
GainSVDTh = 0.001;
isLR = true;
% ----------------------------------------- %

lambda = 1;

% ---- setup psiicos ----- %
SL_rnk = 500;
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

HM = ups.bst.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);

% ---------------------------- proj xyz ------------------------- %
proj_dst_xyz_approx = [0.01, -0.05, 0.1];
proj_dst_ind = ups.FindXYZonGrid(proj_dst_xyz_approx, Ctx_dst.Vertices);
proj_dst_xyz = Ctx_dst.Vertices(proj_dst_ind,:);
proj_dst_indicator = zeros(n_src_dst, 1);
proj_dst_indicator(proj_dst_ind) = 1;
[~, proj_src_ind] = max(Wmat' * proj_dst_indicator);
proj_src_xyz = Ctx_src.Vertices(proj_src_ind,:);
G_proj = HM.gain(:, proj_src_ind * 2 - 1 : proj_src_ind * 2);
% --------------------------------------------------------------- %

% ------------ iflated brain ---------------- %
infl_brain_path = '/home/dmalt/fif_matlab/Brainstorm_db/motorprobe/anat/@default_subject/tess_cortex_pial_high_fig.mat';
Ctx_infl = load(infl_brain_path);

% ------ set seed location on avg brain -------- %
l_M1_dst_xyz = [-0.002, 0.0351, 0.1003];
% l_M1_dst_xyz = [-0.002, -0.0351, 0.1003];
% l_M1_dst_xyz = [0.102, 0.1351, 0.1003];
seed_dst_xyz_approx = l_M1_dst_xyz;
seed_dst_ind = ups.FindXYZonGrid(seed_dst_xyz_approx, Ctx_dst.Vertices);
seed_dst_xyz = Ctx_dst.Vertices(seed_dst_ind,:);
seed_dst_indicator = zeros(n_src_dst, 1);
seed_dst_indicator(seed_dst_ind) = 1;
[~, seed_src_ind]= max(Wmat' * seed_dst_indicator);
seed_src_xyz = Ctx_src.Vertices(seed_src_ind,:);
% ---------------------------------------------- %

threshold = 100;


trials_post = ups.bst.LoadTrials(subj_ID, condition,...
                        freq_band, time_range_post,...
                        HM, protocol_path);
tr_post = trials_post.data;


tr_post_proj = zeros(size(tr_post));
for i_tr = 1:size(tr_post, 3)
    tr_post_proj(:,:,i_tr) = tr_post(:,:,i_tr)...
                           - G_proj * (G_proj' * tr_post(:,:,i_tr));
end

CT_proj = ups.conn.CrossSpectralTimeseries(tr_post_proj, true);


[corr, CT_nosl, Upwr] = ps.PSIICOS(CT_proj, HM.gain, SL_rnk,...
                              sig_rnk, Upwr, seed_src_ind,...
                              cp_part, is_fast);

figure;
h = plot_brain_cmap_hemisplit(Ctx_infl, Ctx_dst, [], corr.data,...
                              zeros(size(corr.data)), 0.1, seed_dst_xyz_approx);
