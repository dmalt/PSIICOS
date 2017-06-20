% ---------------------------------------------------------------- %
% Check if uncoupled sources with high power can generate spurious
% high subspace correlations
% ---------------------------------------------------------------- %
% Date: Sat Jun  3 17:17:35 MSK 2017
% Author: dmalt
% ________________________________________________________________ %

% ------ setup psiicos params ----- %
SL_rnk = 100;
sig_rnk = 20;
Upwr = [];
cp_part = 'real';
is_fast = false;
% --------------------------------- %

% ------- setup simulations ------ %
phase_lag = 0;
n_tr = 100;
GainSVDTh = 0.01;
induced_SNR = 1;
evoked_SNR = 0;
is_use_cache = false;
% --------------------------------- %

[HM, CT, Trials, Ctx, Ctx_HR, XYZGenOut, Ggen] = ups.SimulateData(phase_lag, n_tr, GainSVDTh,...
                                                 induced_SNR, evoked_SNR, is_use_cache);
% --------------------- setup seed ------------------- %
approx_seed_loc =  1.3 * [ 0.05, 0.04, 0.05];
seed_ind = ups.FindXYZonGrid(approx_seed_loc, Ctx.Vertices);
seed_xyz = Ctx.Vertices(seed_ind,:);

n_sen = size(HM.gain, 1);
CT_resh = reshape(mean(CT,2), n_sen, n_sen);
% ------------ remove most powerful sources from data ----------- %
[~, Ps] = ups.conn.DICS(CT_resh, HM.gain, 30, false);

figure;
mask = Ps < 0.9 * max(Ps);
h = plot_brain_cmap_hemisplit(Ctx_HR, Ctx, [], Ps,...
                              zeros(size(Ps)), 0.2, seed_xyz);

filt_ind = reshape([~mask; ~mask]', 1503 * 2, 1);

h = plot_brain_cmap_hemisplit(Ctx_HR, Ctx, [], Ps,...
                              mask, 0.2, seed_xyz);

G_filt = HM.UP * Ggen(:,3:4);

tr_filt = zeros(size(Trials));
n_tr = size(Trials, 3);
for i_tr = 1:n_tr
    tr_filt(:,:, i_tr) = Trials(:,:, i_tr) - G_filt * (pinv(G_filt) * Trials(:,:, i_tr));
end

CT_filt = ups.conn.CrossSpectralTimeseries(tr_filt, false);
CT_filt_resh = reshape(mean(CT_filt,2), n_sen, n_sen);
[~, Ps_filt] = ups.conn.DICS(CT_filt_resh, HM.gain, 30, false);

figure;
mask_filt = zeros(size(Ps_filt));
h = plot_brain_cmap_hemisplit(Ctx_HR, Ctx, [], Ps_filt,...
                              mask_filt, 0.2, seed_xyz);



[corr, CT_post, Upwr] = ps.PSIICOS(CT, HM.gain, SL_rnk,...
                                   sig_rnk, Upwr, seed_ind,...
                                   cp_part, is_fast);
% mask = corr.data < 0.9;
mask_filt = zeros(size(corr.data));
figure;
h = plot_brain_cmap_hemisplit(Ctx_HR, Ctx, [], corr.data,...
                              mask_filt, 0.2, seed_xyz);
