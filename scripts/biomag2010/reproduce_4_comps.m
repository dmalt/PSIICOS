% Bootstrap cross-spetcrum timeseries, then compute SVD
% on real and imaginary parts on each resample and do
% FastScan for first 4 components in real and imag parts
% AUTHOR: dmalt
% DATE Sun Nov 19 22:14:01 MSK 2017
% ______________________________________________________

main
delete(gcp('nocreate'));
current_pool = parpool('local', 4);

time_range = [0.5, 1];
% time_range = [-0.5, 1];

ltr     =(@load_trials);
calc_CT =(@ups.conn.CrossSpectralTimeseries);
proj    =(@ps.ProjectAwayFromPowerComplete);
msvd    =(@svd);

Upwr = ps.GetUpwrComplete(HM.gain, pwr_rnk);
% Upwr = ps.GetUpwrComplete(HM.gain, 101);

IND = ups.indUpperDiag2mat(length(HM.gain)/2);
gpuDevice([]);
con_inds_re = cell(100,4);
% for i_band = 1:length(bands)
for i_band = 1:length(bands)
% freq_band = beta_band;
freq_band = bands{i_band};
% freq_band = alpha_band;

tr_filt= ltr(data_path, time_range, freq_band, HM);
% CT = calc_CT(tr_filt, true);
% CT = ups.GetFakeCT(size(HM.gain,1), 1200);
tic;CTs = ups.bootstrap_CT(tr_filt, 100, true);toc
% profile on;
tic
% i_comp = 1;
parfor i_resamp = 1:length(CTs)
    % tic
    CT_proj = CTs{i_resamp} - Upwr * (Upwr' * CTs{i_resamp});
    % CT_proj = CT - Upwr * (Upwr' * CT);
    [u_re,~,v_re] = msvd(real(CT_proj));
    [u_im,~,v_im] = msvd(imag(CT_proj));
    % [u,s,v] = msvd(imag(CT_proj));
    % [u,s,v] = msvd(imag(CT));
    % profile on
    for i_comp = 1:4
        vs_re{i_resamp}{i_comp} = v_re(:,i_comp);
        vs_im{i_resamp}{i_comp} = v_im(:,i_comp);
        % [CS_re, ~] = ps.PSIICOS_ScanFast(HM.gain, u_re(:,i_comp));
        % [CS_im, IND] = ps.PSIICOS_ScanFast(HM.gain, u_im(:,i_comp));
        % gd = gpuDevice;
        % idx = gd.Index;
        % disp(['Using GPU', num2str(idx)]);
        % tic
        CS_re = PSIICOS_ScanFastGPU(HM.gain, u_re(:,i_comp), false, 2000);
        con_inds_re{i_resamp}{i_comp} = ups.threshold_connections(CS_re, threshold, IND);
        CS_re = [];
        CS_im = PSIICOS_ScanFastGPU(HM.gain, u_im(:,i_comp), false, 2000);
        con_inds_im{i_resamp}{i_comp} = ups.threshold_connections(CS_im, threshold, IND);
        CS_im = [];
        % toc
        % con_inds_re{i_comp}{i_resamp} = ups.threshold_connections(CS_re, threshold, IND);
    end
    % profile viewer
    % disp(['resample number is ', num2str(i_resamp)]);
    % toc
end
toc
save([band_names{i_band}, '.mat'], 'vs_re', 'vs_im', 'con_inds_re', 'con_inds_im', 'CTs', '-v7.3');
end

delete(current_pool)
return;
%  Things to check bootstrap solution {{{1 %
    con = ups.Bundles(con_inds_re{1}, HM, CtxInfl);
    con.Plot()
    con_av = con.Average();
    i = con_av.conInds{1}(1);
    j = con_av.conInds{1}(2);
    gi = con_av.headModel.gain(:, i * 2 - 1 : 2 * i);
    gj = con_av.headModel.gain(:, j * 2 - 1 : 2 * j);
    Gij = kron(gi,gj);
    Gji = kron(gj,gi);
    GG = [Gji, Gij];
    [u,s,v] = svd(GG);
    rnk_range =  1:8;
    CT_proj_2 = CT_proj_1 - u(:,rnk_range) * (u(:,rnk_range)' * CT_proj_1);

    [u_re,~,~] = msvd(real(CT_proj_2));
    [CS_re_2, IND] = ps.PSIICOS_ScanFast(HM.gain, u_re(:,i_comp));

    con_inds_re{i_comp} = ups.threshold_connections(CS_re_2, threshold, IND);


    threshold = 30;
    ch_path = [protocol_path, '/data/biomag2010/@rawdataset02/channel_ctf_acc1.mat'];
    ch = load(ch_path);

    i_band = 3;
    freq_band = bands{i_band};
    tr_filt= ltr(data_path, time_range, freq_band, HM);
    CT = calc_CT(tr_filt, true);

    i_comp = 1;

    rn = 100;
    Upwr = ps.GetUpwrComplete(HM.gain, rn);
    CT_proj = CT - Upwr * (Upwr' * CT);
    [u_re,s_re,v_re] = svd(real(CT_proj));
    % [u_re,s_re,v_re] = svd(real(CT));
    [CS_re, IND] = ps.PSIICOS_ScanFast(HM.gain, u_re(:,i_comp));
    con_inds_re{i_comp} = ups.threshold_connections(CS_re, threshold, IND);
    con = ups.Bundles(con_inds_re{i_comp}, HM, CtxInfl);
    figure('Name', num2str(rn));
    con.Plot();

    uu = u_re(:,i_comp);
    uu = reshape(uu,33,33);
    uu = HM.UP' * uu * HM.UP;
    uu = abs(uu);
    [i,j] = find(uu>0.02);
    ups.plt.DrawConnectionsOnSensors([i,j], ch_path, true);

%     con_m = con.Merge();
%     cc = con_m.Clusterize(1, 0.0001);
%     cca  = cc.Average();
%     ccam = cca.Merge();
%     ccaa = ccam.Average();
%     cca.conInds = {ccaa.conInds{:}, cca.conInds{:}};

%     n_clust = length(cc.conInds);
%     sizes = zeros(n_clust,1);
%     linewidth = zeros(n_clust,1);

%     for i_clust = 1:length(cca.conInds)
%         sizes(i_clust) = size(cc.conInds{i_clust},1);
%     end



%     linewidth =  floor(sizes / 30) + 1;
%     linewidth = ones(size(sizes));
%     linewidth_ = [12; linewidth]
%     m_radii = 0.001 * (1 + 0.4 * log(sizes));
%     m_radii = ones(size(sizes)) * 0.001;
%     m_radii_ = [0.005; m_radii];
%     icols = ones(size(m_radii_));
%     icols(1) = 2;
%     alpha = 0.2 * ones(size(m_radii_));
%     alpha(1) = 1;
%     cca.PlotViews(0.2, linewidth_, m_radii_, icols, alpha)
%     % % con.CT = CT_proj;
%     % con.CT = CT_proj;
%     % con_av = con.Average();
%     % con_av.PlotViews(0.2);


%     % vv = [vs{:}];
%     % m_v = mean(vv,2);
%     % std_v = std(vv,[],2) / 10;
%     % errorbar(m_v, std_v);
%     % figure; con.PlotViews(0.1, 2);


%     % plot(v(1,:))
%     % [corr, Cpvec, Upwr] =
%     % colors = ups.plt.GetColors(3);

%     % figure;
%     % c = zeros(4,3);
%     % for i = 1:4
%     %     c(i,:) = colors(i + 1).RGB;
%     %     plot(v(:,i), 'LineWidth', 2, 'Color', colors(i + 1).RGB);
%     %     hold on;
%     % end

%     % legend('comp1', 'comp2', 'comp3', 'comp4');
%  1}}} %
