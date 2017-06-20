% ---------------------------------------------------------------- %
% Check if uncoupled sources with high power can generate spurious
% high subspace correlations
% ---------------------------------------------------------------- %
% Date: Sat Jun  3 17:17:35 MSK 2017
% Author: dmalt
% ________________________________________________________________ %

% ------ setup psiicos params ----- %
SL_rnk = 350;
sig_rnk = 20;
Upwr = [];
cp_part = 'real';
is_fast = true;
% --------------------------------- %

% ------- setup simulations ------ %
phase_lag = 0;
n_tr = 100;
GainSVDTh = 0.01;
induced_SNR = 10;
evoked_SNR = 0;
is_use_cache = false;
% --------------------------------- %

[HM, CT, Trials, Ctx, Ctx_HR] = ups.SimulateData(phase_lag, n_tr, GainSVDTh,...
                                                 induced_SNR, evoked_SNR, is_use_cache);

% ----------- setup seed --------- %
approx_seed_loc =  1.3 * [ 0.05,  0.04, 0.05];
seed_ind = ups.FindXYZonGrid(approx_seed_loc, Ctx.Vertices);
seed_xyz = Ctx.Vertices(seed_ind,:);


[corr, CT_post, Upwr] = ps.PSIICOS(CT, HM.gain, SL_rnk,...
                                sig_rnk, Upwr, seed_ind,...
                                cp_part, is_fast);
mask = corr.data < 0.9;
figure;
h = plot_brain_cmap_hemisplit(Ctx_HR, Ctx, [], corr.data,...
                              mask, 0.2, seed_xyz);
