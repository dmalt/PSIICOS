% Preparation steps for analyses
% AUTHOR: dmalt
% ------------------------------ %

% -------- setup  constants ------- %
% paths
home = getenv('HOME');
bst_path = [home,'/Documents/MATLAB/bst/brainstorm_db'];
data_path = [home, '/Data/mentrot/MentalRotationDeLange/preprocessed'];
subj_ID = 'biomag2010';
protocol_path = [bst_path,'/mentrot'];
condition = 'raw';

% head model
isLR = false;
% isLR = false;
GainSVDTh = 0.001; % results in 45 components
ch_type = 'MEG';


% load head model
HM = ups.bst.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh, ch_type);



% ch_path = '/home/dmalt/bst/brainstorm_db/mentrot/data/biomag2010/@rawdataset02/channel_ctf_acc1.mat';
% ch = load(ch_path);
% ch = ch.Channel;
% ch_used = ups.bst.PickChannels(ch, 'MEG');
% ch_meg = ch(ch_used);
% ups.plt.DisplayTopographyCTF(pwr, ch_meg);

% lambda = 1;
% [A, Ps] = ups.conn.DICS(CT_resh, HM.gain, lambda);
[Ctx, CtxHR, CtxHHR] = ups.bst.GetCtx(subj_ID, protocol_path);
% figure;
% plot_brain_cmap_hemisplit(CtxHR, Ctx, [], Ps, zeros(size(Ps)), 0.2);
% % plot_brain_cmap(CtxHHR, Ctx, [], Ps, zeros(size(Ps)), 0.2);

pwr_rnk = 150;
threshold = 30;
% CT_proj_from_vc = ps.ProjectAwayFromPowerFixedOr(CT_resh(:), HM.gain, pwr_rnk);
% [CS, IND] = ps.PSIICOS_ScanFast(HM.gain, real(CT_proj_from_vc));
% con_inds = ups.threshold_connections(CS, threshold, IND);

% con = ups.Bundles(con_inds, HM, CtxHHR);
% % con_c = con.Clusterize(10,0.02);
% figure;
% con.Plot(0.2);
