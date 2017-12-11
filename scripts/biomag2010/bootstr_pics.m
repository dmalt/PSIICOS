% Draw pictures for real data bootstraps
% -------------------------------------- %
% AUTHOR: dmalt
% DATE: Thu Nov 16 13:48:51 MSK 2017
% ______________________________________ %

main;
data_dir = './bootstrap_data/';
save_dir = './pics/bootstrap_pics/';
n_comps = 4;
% i_band = 3;
% i_comp = 2;
% for i_band = 1:length(band_names)
for i_band = 7:length(band_names)
% for i_band = 3:3
% for i_band = 4:4%length(band_names)
    s = load([data_dir, band_names{i_band}, '.mat']);
    % t = load([data_dir, band_names{i_band}, '_fix.mat']);
    for i_comp = 1:n_comps
    % for i_comp = 1:2
        con_re = ups.Bundles(s.con_inds_re{i_comp}, HM, CtxInfl);
        % con_re = ups.Bundles(t.con_inds_re{i_comp}, HM, CtxInfl);
        % cc = ups.Bundles(s.con_inds_re{i_comp}, HM, CtxHHR);
        con_re_av = con_re.Average(false);
        con_re_av_clust = con_re_av.Clusterize(0.01);
        % con_re_big = con_re_av_clust.Threshold(10);
        % con_re_centers = con_re_big.Average();
        cc = con_re_av_clust([]);

        for i_cl = 1:length(con_re_av_clust.conInds)
            % ---- merge superimposed sticks ----- %
            tt = con_re_av_clust(i_cl).Clusterize(0.0001);
            ttt = GravitateAlpha(tt);
            ttt.icol(:) = con_re_av_clust.icol(i_cl);
            cc = cc + ttt;
        end

        con_re_final = cc;
        figname = [band_names{i_band}, '_', num2str(i_comp), '_', 're'];
            % figure('Name', figname);
        con_re_final.PlotViews(0.2);
        set(gcf, 'Color','None');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        % pause(2);
        export_fig([save_dir, figname, '.png'], '-png', '-transparent');
        close(gcf);

        con_im = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxInfl);
        % cc = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxHHR);
        con_im_av = con_im.Average(false);
        con_im_av_clust = con_im_av.Clusterize(0.01);
        % con_im_big = con_im_av_clust.Threshold(10);
        % con_im_centers = con_im_big.Average();
        cc = con_re_av_clust([]);
        for i_cl = 1:length(con_im_av_clust.conInds)
            % ---- merge superimposed sticks ----- %
            tt = con_im_av_clust(i_cl).Clusterize(0.0001);
            ttt = GravitateAlpha(tt);
            ttt.icol(:) = con_im_av_clust.icol(i_cl);
            cc = cc + ttt;
        end

        % con_im_centers.m_rad(:) = 0.004;
        % con_im_centers.lwidth(:) = 4;
        % con_im_av_clust.alpha(:) = 0.4;
        % con_im_av_clust.lwidth(:) = 1;
        % con_im_av_clust.m_rad(:) = 0.0015;
        % con_im_final = con_im_av_clust + con_im_centers;
        con_im_final = cc;
        figname = [band_names{i_band}, '_', num2str(i_comp), '_', 'im'];
            % figure('Name', figname);
        con_im_final.PlotViews(0.2);
        set(gcf, 'Color','None');
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        % pause(2);
        export_fig([save_dir, figname, '.png'], '-png', '-transparent');
        close(gcf);
        % con_im = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxHHR);
        % % cc = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxHHR);
        % con_im_av = con_im.Average(false);
        % con_im_av_clust = con_im_av.Clusterize(0.01);
        % con_im_big = con_im_av_clust.Threshold(10);
        % con_im_centers = con_im_big.Average();
        % % con_im_centers.m_rad(:) = 0.004;
        % % con_im_centers.lwidth(:) = 4;
        % con_im_av_clust.alpha(:) = 0.4;
        % con_im_av_clust.lwidth(:) = 1;
        % con_im_av_clust.m_rad(:) = 0.0015;
        % con_im_final = con_im_av_clust + con_im_centers;
        % figname = [band_names{i_band}, '_', num2str(i_comp), '_', 'im'];
        % % figure('Name', figname);
        % con_im_final.PlotViews(0.2);
        % set(gcf, 'Color','None');
        % % pause(2);
        % export_fig([save_dir, figname, '.png'], '-png', '-transparent');
        % close(gcf);
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
end

