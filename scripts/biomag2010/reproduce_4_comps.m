main

time_range = [0, 1];

alpha_band = [8,12];
beta_band = [16,24];
gamma_band = [65, 85];
lowgamma_band = [30, 60];
theta_band = [4,8];
bands = {theta_band, alpha_band, beta_band, lowgamma_band, gamma_band};
band_names = {'theta_band', 'alpha_band', 'beta_band', 'lowgamma_band', 'gamma_band'};


ltr     = memoize(@load_trials);
calc_CT = memoize(@ups.conn.CrossSpectralTimeseries);
proj    = memoize(@ps.ProjectAwayFromPowerComplete);
msvd    = memoize(@svd);

Upwr = ps.GetUpwrComplete(HM.gain, pwr_rnk);

for i_band = 1:length(bands)
% freq_band = beta_band;
freq_band = bands{i_band};
% freq_band = alpha_band;

tr_filt= ltr(data_path, time_range, freq_band, HM);
% CT = calc_CT(tr_filt(:,:,10:15), true);
% CT = ups.GetFakeCT(size(HM.gain,1), 1200);
CTs = ups.bootstrap_CT(tr_filt, 100, true);
% profile on;
for i_comp = 1:4
% i_comp = 1;
    for i_resamp = 1:length(CTs)
        tic
        CT_proj = CTs{i_resamp} - Upwr * (Upwr' * CTs{i_resamp});
        [u_re,s_re,v_re] = msvd(real(CT_proj));
        [u_im,s_im,v_im] = msvd(imag(CT_proj));
        % [u,s,v] = msvd(imag(CT_proj));
        % [u,s,v] = msvd(imag(CT));
        vs_re{i_comp}{i_resamp} = v_re(:,i_comp);
        vs_im{i_comp}{i_resamp} = v_im(:,i_comp);
        [CS_re, IND] = ps.PSIICOS_ScanFast(HM.gain, u_re(:,i_comp));
        [CS_im, IND] = ps.PSIICOS_ScanFast(HM.gain, u_im(:,i_comp));
        con_inds_re{i_comp}{i_resamp} = ups.threshold_connections(CS_re, threshold, IND);
        con_inds_im{i_comp}{i_resamp} = ups.threshold_connections(CS_im, threshold, IND);
        toc
    end
end
save([band_names{i_band}, '.mat'], 'vs_re', 'vs_im', 'con_inds_re', 'con_inds_im', '-v7.3');
end

%     con = ups.Connections(con_inds, HM, CtxHHR);
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
