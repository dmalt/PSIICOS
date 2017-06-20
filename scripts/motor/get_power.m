lambda = 1;
G = HM.gain;
n_sen = size(G,1);
C_pre = reshape(mean(CT_pre, 2), n_sen, n_sen);
C_post = reshape(mean(CT_post, 2), n_sen, n_sen);


[A, Ps_pre, Cs, IND] = ups.conn.DICS(C_pre, G, lambda);
[A, Ps_post, Cs, IND] = ups.conn.DICS(C_post, G, lambda);


figure;
h = plot_brain_cmap_hemisplit(CtxHHR_dst, Ctx_dst, [], Wmat * (Ps_pre -  Ps_post),...
                              zeros(size(Wmat * Ps)), 0.1, seed_dst_xyz);
