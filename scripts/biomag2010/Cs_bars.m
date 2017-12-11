% Plot Cs values from ScanFast for each frequency band
% and for each SVD component as bars to see if the higher values
% of Cs on individual runs correspond to better stability on bootstrap
%
% Summary:
%   No monotonicity in max(Cs) dependence on the number of svd component
%   whatsoever. Real part max(Cs) shows one order of magnitude smaller
%   values than that for the imaginary part.
% --------------------------------------------------------------------
% AUTHOR: dmalt
% DATE: Sun Nov 19 23:29:25 MSK 2017
% ____________________________________________________________________
main;

% n_comps = 4;

bands = {theta_band, alpha_band, beta_band, lowgamma_band, gamma_band};
band_names = {'theta_band', 'alpha_band', 'beta_band', 'lowgamma_band', 'gamma_band'};
ltr     = memoize(@load_trials);
calc_CT = memoize(@ups.conn.CrossSpectralTimeseries);
% proj    = memoize(@ps.ProjectAwayFromPowerComplete);
msvd    = memoize(@svd);

% Upwr = ps.GetUpwrComplete(HM.gain, pwr_rnk);

% for i_band = 2:length(band_names)
    % freq_band = bands{i_band};
    % tr_filt= ltr(data_path, time_range, freq_band, HM);
    % CT = calc_CT(tr_filt, true);
%     for i_comp = 1:n_comps
%         tic
        % CT_proj = CT - Upwr * (Upwr' * CT);
        % [u_re,s_re,v_re] = msvd(real(CT_proj));
        % [u_im,s_im,v_im] = msvd(imag(CT_proj));
%         vs_re{i_comp} = v_re(:,i_comp);
%         vs_im{i_comp} = v_im(:,i_comp);
%         [CS_re{i_comp}, IND] = ps.PSIICOS_ScanFast(HM.gain, u_re(:,i_comp));
%         [CS_im{i_comp}, IND] = ps.PSIICOS_ScanFast(HM.gain, u_im(:,i_comp));
%         con_inds_re{i_comp} = ups.threshold_connections(CS_re{i_comp}, threshold, IND);
%         con_inds_im{i_comp} = ups.threshold_connections(CS_im{i_comp}, threshold, IND);
%         toc
%     end
% save([band_names{i_band}, '.mat'], 'vs_re', 'vs_im', 'con_inds_re', 'con_inds_im', 'CS_re', 'CS_im', '-v7.3');
% end
% ij = con_re_centers.conInds{1};
% gi = HM.gain(:, ij(1) * 2 - 1 : ij(1) * 2);
% gj = HM.gain(:, ij(2) * 2 - 1 : ij(2) * 2);
% Gij = kron(gi,gj);


main
data_dir = './Cs_data/';
save_dir = './pics/Cs_pics/';
n_comps = 4;
i_band = 1;

% for i_band = 1:length(band_names)
% for i_band = 4:4%length(band_names)
    s = load([data_dir, band_names{i_band}, '.mat']);
    for i_comp = 1:n_comps
        con_re = ups.Bundles(s.con_inds_re{i_comp}, HM, CtxHHR);
        % cc = ups.Bundles(s.con_inds_re{i_comp}, HM, CtxHHR);
        con_re_clust = con_re.Clusterize(0.01);
        con_re_big = con_re_clust.Threshold(5);
        con_re_centers = con_re_big.Average();
        % con_re_centers.m_rad(:) = 0.004;
        % con_re_centers.lwidth(:) = 4;
        con_re_av_clust.alpha(:) = 0.4;
        con_re_av_clust.lwidth(:) = 1;
        con_re_av_clust.m_rad(:) = 0.001;
        con_re_final = con_re_av_clust + con_re_centers;
        figname = [band_names{i_band}, '_', num2str(i_comp), '_', 're'];
        % figure('Name', figname);
        con_re_final.PlotViews(0.2);
        set(gcf, 'Color','None');
        % pause(2);
        export_fig([save_dir, figname, '.png'], '-png', '-transparent');
        close(gcf);

        con_im = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxHHR);
        % cc = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxHHR);
        con_im_av = con_im.Average(false);
        con_im_av_clust = con_im_av.Clusterize(0.01);
        con_im_big = con_im_av_clust.Threshold(10);
        con_im_centers = con_im_big.Average();
        % con_im_centers.m_rad(:) = 0.004;
        % con_im_centers.lwidth(:) = 4;
        con_im_av_clust.alpha(:) = 0.4;
        con_im_av_clust.lwidth(:) = 1;
        con_im_av_clust.m_rad(:) = 0.001;
        con_im_final = con_im_av_clust + con_im_centers;
        figname = [band_names{i_band}, '_', num2str(i_comp), '_', 'im'];
        % figure('Name', figname);
        con_im_final.PlotViews(0.2);
        set(gcf, 'Color','None');
        % pause(2);
        export_fig([save_dir, figname, '.png'], '-png', '-transparent');
        close(gcf);
        % con_im = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxHHR);
        % con_im_av = con_im.Average();
        % con_im_av_clust = con_im_av.Clusterize(1, 0.01);
        % figname = [band_names{i_band}, '_', num2str(i_comp), '_', 'im'];
        % figure('Name', figname);
        % con_im_av_clust.PlotViews(0.1, 1, 0.001, []);
        % pause(2);
        % set(gcf, 'Color','None');
        % export_fig([save_dir, figname, '.png'], '-png', '-transparent');
        % close(gcf);
    end
% end


% n_svd = size(u_re,2);
% scorr = zeros(n_svd,1);
% for i = 1:n_svd
%     scorr(i) = ps.subcorr([Gij, Gji], u_im(:,i));
% end
