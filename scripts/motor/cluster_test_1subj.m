
protocol_path = '/home/dmalt/fif_matlab/Brainstorm_db/motorprobe';
subj_ID = 'AB';
suffix = '_control_RH_raw_tsss_ica';

condition = [subj_ID, suffix];
freq_band = [19, 25];
time_range_pre = [-0.6, -0.2];
time_range_post = [0, 0.4];

n_resamp = 5000;

% ---- setup head model loading params ---- %
GainSVDTh = 0.01;
isLR = true;
% ----------------------------------------- %

lambda = 1;

% ---- setup psiicos ----- %
SL_rnk = 350;
sig_rnk = 0;
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


HM = ups.bst.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);

[tr_pre, ~, Ps_pre_dst,...
 CS_pre_dst, Upwr, Wmat] = prepare_cond(subj_ID, condition,...
                                        protocol_path,...
                                        HM, freq_band,...
                                        time_range_pre,...
                                        lambda, seed_src_ind,...
                                        SL_rnk, sig_rnk,...
                                        Upwr, cp_part,...
                                        is_fast, Wmat);

[tr_post, ~, Ps_post_dst, CS_post_dst, ~, ~] = prepare_cond(subj_ID, condition,...
                                                      protocol_path,...
                                                      HM, freq_band,...
                                                      time_range_post, lambda,...
                                                      seed_src_ind, SL_rnk, sig_rnk,...
                                                      Upwr, cp_part, is_fast, Wmat);
CS_dst = (CS_post_dst.data - CS_pre_dst.data);


clust_threshold = 0.12;

CS_dst_perm = zeros(n_src_dst, n_resamp);
adj_mat = get_adjacency_matrix(Ctx_dst);
clust_stat = cell(n_resamp,1);
for i_resamp = 1:n_resamp
    tic;
    [tr_part1, tr_part2] = get_trials_partition(tr_pre, tr_post);
    CT_part1 = ups.conn.CrossSpectralTimeseries(tr_part1, true);
    CT_part2 = ups.conn.CrossSpectralTimeseries(tr_part2, true);
    corr1 = ps.PSIICOS(CT_part1, HM.gain, SL_rnk, sig_rnk,...
                       Upwr, seed_src_ind, cp_part, is_fast);
    CS1_dst_perm = Wmat * corr1.data;

    corr2 = ps.PSIICOS(CT_part2, HM.gain, SL_rnk, sig_rnk,...
                    Upwr, seed_src_ind, cp_part, is_fast);
    CS2_dst_perm = Wmat * corr2.data;

    % difference null-distribution
    CS_dst_perm(:,i_resamp) = CS2_dst_perm - CS1_dst_perm;
    clust_mask_inc = CS_dst_perm(:,i_resamp) > clust_threshold;
    clust_mask_dec = CS_dst_perm(:,i_resamp) < -clust_threshold;

    [clusters_inc, clusters_dec] = get_clusters(adj_mat, clust_mask_inc, clust_mask_dec);


    stat_inc = zeros(length(clusters_inc),1);
    for iclust = 1:length(clusters_inc)
        stat_inc(iclust) = sum(CS_dst_perm(clusters_inc{iclust}, i_resamp));
    end

    stat_dec = zeros(length(clusters_dec),1);
    for iclust = 1:length(clusters_dec)
        stat_dec(iclust) = abs(sum(CS_dst_perm(clusters_dec{iclust}, i_resamp)));
    end


    clust_stat{i_resamp} = sort([stat_inc; stat_dec],1,'descend');
end

% ----------- get preliminary p-values -------------- %
p_temp = ones(n_src_dst,1);
for i_src = 1:n_src_dst
    p_temp(i_src) = (sum(CS_dst_perm(i_src,:) > abs(CS_dst(i_src)))...
                  + sum(CS_dst_perm(i_src,:) < -abs(CS_dst(i_src))))...
                  / n_resamp;
end

stat_mask = p_temp > 0.03;



% CS_dst = (CS_post_dst - CS_pre_dst) ./ CS_pre_dst;

figure;
h = plot_brain_cmap_hemisplit(CtxHHR_dst, Ctx_dst, [], p_temp,...
                              stat_mask, 0.1, seed_dst_xyz_approx);

