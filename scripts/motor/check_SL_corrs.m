% NEED:
% seed_src_ind
% CT_post
% HM
% Upwr
CT_proj = CT_post;
[u,s,v] = svd(real(CT_proj));
UU = u(:,20);

G = HM.gain;
[n_sen, n_src] = size(G);
n_src = n_src / 2;

A = zeros(n_sen ^ 2, 3);
gg_seed = zeros(n_sen ^ 2, 3);

g_seed_i = G(:, seed_src_ind * 2 - 1);
g_seed_j = G(:, seed_src_ind * 2);

gg_seed(:,1) = kron(g_seed_i, g_seed_i);
gg_seed(:,2) = kron(g_seed_j, g_seed_j);
gg_seed(:,3) = kron(g_seed_i, g_seed_j) + kron(g_seed_j, g_seed_i);

gg_seed = gg_seed - Upwr * (Upwr' * gg_seed);

corr = zeros(n_src,1);

for i_src = 1:n_src
    g_i = G(:, i_src * 2 - 1);
    g_j = G(:, i_src * 2);

    A(:, 1) = kron(g_i, g_i);
    A(:, 2) = kron(g_j, g_j);
    A(:, 3) = kron(g_i, g_j) + kron(g_j, g_i);
    A = A - Upwr * (Upwr' * A);

    sl_subs = [gg_seed, A];
    c = ps.subcorr(sl_subs, UU);
    corr(i_src) = c(1);
end

figure;
h = plot_brain_cmap_hemisplit(CtxHHR_dst, Ctx_dst, [], corr,...
                              zeros(size(corr)), 0.1, seed_dst_xyz_approx);


