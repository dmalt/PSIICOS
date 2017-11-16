% Draw pictures for real data bootstraps
% -------------------------------------- %
% AUTHOR: dmalt
% DATE: Thu Nov 16 13:48:51 MSK 2017
% ______________________________________ %

main;
band_names = {'theta_band', 'alpha_band', 'beta_band', 'lowgamma_band', 'gamma_band'};
data_dir = './bootstrap_data/';
save_dir = './pics/bootstrap_pics/';
n_comps = 4;

for i_band = 1:length(band_names)
% for i_band = 4:4%length(band_names)
    s = load([data_dir, band_names{i_band}, '.mat']);
    for i_comp = 1:n_comps
        con_re = ups.Connections(s.con_inds_re{i_comp}, HM, CtxHHR);
        con_re_av = con_re.Average();
        con_re_av_clust = con_re_av.Clusterize(1, 0.01);
        figname = [band_names{i_band}, '_', num2str(i_comp), '_', 're'];
        figure('Name', figname);
        con_re_av_clust.PlotViews(0.1, 1, 0.001, []);
        set(gcf, 'Color','None');
        pause(2);
        export_fig([save_dir, figname, '.png'], '-png', '-transparent');
        close(gcf);

        con_im = ups.Connections(s.con_inds_im{i_comp}, HM, CtxHHR);
        con_im_av = con_im.Average();
        con_im_av_clust = con_im_av.Clusterize(1, 0.01);
        figname = [band_names{i_band}, '_', num2str(i_comp), '_', 'im'];
        figure('Name', figname);
        con_im_av_clust.PlotViews(0.1, 1, 0.001, []);
        pause(2);
        set(gcf, 'Color','None');
        export_fig([save_dir, figname, '.png'], '-png', '-transparent');
        close(gcf);
    end
end

