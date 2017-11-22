% Plot beard pics
% -------------------------------------- %
% AUTHOR: dmalt
% ______________________________________ %

% main;
band_names = {'theta_band', 'alpha_band', 'beta_band', 'lowgamma_band', 'gamma_band'};
data_dir = './bootstrap_data/';
save_dir = './pics/bootstrap_pics/';
n_comps = 4;
i_band = 3;
i_comp = 1;

% for i_band = 1:length(band_names)
% for i_band = 4:4%length(band_names)
    s = load([data_dir, band_names{i_band}, '.mat']);
    for i_comp = 1:n_comps
        v_re = [s.vs_re{i_comp}{:};];
        v_im = [s.vs_im{i_comp}{:};];

        v_re_mean = mean(v_re,2);
        v_im_mean = mean(v_im,2);

        v_re_std = std(v_re,[],2) / 10;
        v_im_std = std(v_im,[],2) / 10;

        errorbar(v_re_mean(20:end-20), v_re_std(20:end-20))
        % plot(v_re_mean(20:end-20))
        hold on;
    end
% end

