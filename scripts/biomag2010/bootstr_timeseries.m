% Estimate timeseries for bootstrap clusters
% -------------------------------------- %
% AUTHOR: dmalt
% DATE: Sat Dec 2 16:58:08 MSK 2017
% ______________________________________ %

main;

data_dir = './bootstrap_data/';
save_dir = './pics/bootstrap_pics/';
n_comps = 4;
i_band = 1;
i_comp = 1;
for i_band = 7:length(band_names)
    s = load([data_dir, band_names{i_band}, '.mat']);
    % t = load([data_dir, band_names{i_band}, '_fix.mat']);
    for i_comp = 1:n_comps
        % con_re = ups.Bundles(t.con_inds_re{i_comp}, HM, CtxInfl);
        con_re = ups.Bundles(s.con_inds_re{i_comp}, HM, CtxInfl);
        con_re_av = con_re.Average(false);
        con_re_av_clust = con_re_av.Clusterize(0.01);
        plot_clust_tseries(con_re_av_clust, con_re_av, HM,...
                           s.CTs, band_names{i_band}, 'Re', i_comp);
% return

        con_im = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxInfl);
        % cc = ups.Bundles(s.con_inds_im{i_comp}, HM, CtxHHR);
        con_im_av = con_im.Average(false);
        con_im_av_clust = con_im_av.Clusterize(0.01);
        plot_clust_tseries(con_im_av_clust, con_im_av, HM,...
                           s.CTs, band_names{i_band}, 'Im', i_comp);
    end
end

% 3630 9698
