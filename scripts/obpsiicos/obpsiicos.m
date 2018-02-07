%  setup input {{{1 %
subjID = 'test';
PhaseLag = pi / 20;
% PhaseLag = pi / 2 - pi / 20;
GainSVDTh = 0.01;
sim_cached = memoize(@ups.SimulateData);
[HM, CT, Trials, Ctx, XYZGen] = sim_cached(PhaseLag, 100, GainSVDTh, 0.3, 0, false);
G2d = HM.gain;
n_ch = size(G2d, 1);
r_g = rank(G2d);
rr = r_g * (r_g - 1) / 2 + r_g;

tril_mask = true(r_g);
tril_mask = tril(tril_mask);
tril_mask = tril_mask(:);

% fname = 'headmodel_vol_os_meg.mat';
% gg = load(fname);
% G = gg.Gain(36:186,:);
% G2 = ups.ReduceToTangentSpace(G);
% G2d = ups.ReduceSensorSpace(G2, GainSVDTh);

%  1}}} %


% profile on;
% for i = 1:n_sites
tic
C_sl = compute_re_cov_tril(G2d);
C_sl = double(C_sl);
W = sqrtm(C_sl);
W_inv = inv(W);
toc
% profile viewer

CT_tril_re = real(CT(tril_mask,:));

% %  compute whitener {{{2 %
% [u, s, v] = svd(C_sl);
% ss = diag(s);
% ssw_re = ones(rr, 1) ./ sqrt(ss(1 : rr));
% ssw_im = zeros(r_g ^ 2 - rr, 1);
% ssw = [ssw_re; ssw_im];
% W_inv = u * diag(ssw) * v';
% W = u * diag([sqrt(ss(1:rr)); ssw_im]) * v';

% lambda = 10;
% C_sl_reg = (C_sl + lambda * trace(C_sl) / n_ch ^ 2 * eye(size(C_sl)));
% W = sqrtm(C_sl_reg);

tic;
% P = obProjector(G2d, 350, W_inv);
P = proj_tril(G2d, 400, W_inv);
toc;
%  2}}} %

% W = sqrtm(C_sl);

% %  whiten data and find sources{{{3 %
% % CT_w = inv(W) *  CT;
% CT_w = W_inv *  CT;
% CT_pw = CT_w - P * (P' * CT_w);
% CT_wpw = W * CT_pw;
% [Cs_obp, IND] = ps.PSIICOS_ScanFast(G2d, mean(CT_wpw,2));

CT_w = W_inv *  CT_tril_re;
CT_pw = CT_w - P * (P' * CT_w);
CT_wpw = W * CT_pw;
CT_rest = zeros(size(CT));
CT_rest(tril_mask,:) = CT_wpw;
for i_time = 1:size(CT_rest, 2)
    CC = reshape(CT_rest(:,i_time), r_g, r_g);
    CC = CC + CC' - diag(diag(CC));
    CT_rest(:,i_time) = CC(:);
end
CT_rest = CT_rest + 1i * imag(CT);
[Cs_obp, IND] = ps.PSIICOS_ScanFast(G2d, (mean(CT_rest,2)));
ij_ps = ups.threshold_connections(Cs_obp, 200, IND);
con_ps = ups.Bundles(ij_ps, HM, Ctx);
figure;
con_ps.Plot();
%  3}}} %

%  look at original psiicos {{{4 %
CT_p = ps.ProjectAwayFromPowerComplete(CT, G2d, 200);
[Cs_p, IND] = ps.PSIICOS_ScanFast(G2d, mean(CT_p, 2));
ij_ps = ups.threshold_connections(Cs_p, 200, IND);
con_ps = ups.Bundles(ij_ps, HM, Ctx);
figure;
con_ps.Plot();
%  4}}} %

%  plot connections {{{5 %
%  5}}} %

% con_ps_clust = con_ps.Clusterize(10, 0.02);
% con_ps_av = con_ps_clust.Average();
% con_ps_av.Plot(0.2, 4, 0.004);


%  generate ROC curve {{{6 %
NetworkPairIndex{1} = [1,2];
NetworkPairIndex{2} = [1,2,3];
NPI = NetworkPairIndex{2};
Dmax = 0.02;
n_levels = 200;
GridLoc = HM.GridLoc;
[SPC_obp, TPR_obp, PPV_obp] = ups.GenerateScores(Cs_obp, Dmax, GridLoc, IND, n_levels, XYZGen, NPI);
[SPC_p, TPR_p, PPV_p] = ups.GenerateScores(Cs_p, Dmax, GridLoc, IND, n_levels, XYZGen, NPI);
figure;
plot(1 - SPC_obp, TPR_obp);
hold on;
plot(1 - SPC_p, TPR_p);
legend('OBPSIICOS', 'PSIICOS')
figure;
plot(TPR_obp, PPV_obp);
hold on;
plot(TPR_p, PPV_p);
legend('OBPSIICOS', 'PSIICOS');
%  6}}} %
