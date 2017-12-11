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

n_times = 1200;
n_comps = 4;

v_re_mean = zeros(n_times,n_comps);
v_im_mean = zeros(n_times,n_comps);
v_re_std = zeros(n_times,n_comps);
v_im_std = zeros(n_times,n_comps);


% -------------- define ticks ------------- %
sfreq = 1200;
ticks = 0:150:1200;
label_vals = ticks / sfreq * 1000;
label_names = cell(size(label_vals));
for ii = 1:length(label_vals)
    label_names{ii} = num2str(label_vals(ii));
end
% ----------------------------------------- %

for i_band = 1:length(band_names)
% for i_band = 4:4%length(band_names)
    s = load([data_dir, band_names{i_band}, '.mat']);
    for i_comp = 1:n_comps
        v_re = [s.vs_re{i_comp}{:};];
        v_im = [s.vs_im{i_comp}{:};];

        v_re = flip_components(v_re);
        v_im = flip_components(v_im);

        v_re_mean(:,i_comp) = mean(v_re,2);
        v_im_mean(:,i_comp) = mean(v_im,2);

        v_re_std(:,i_comp) = std(v_re,[],2) / 10;
        v_im_std(:,i_comp) = std(v_im,[],2) / 10;

        % plot(v_re_mean(20:end-20))
        % hold on;
    end
    figname = ['tseries_', band_names{i_band}];
    figure('Name', figname);
    % subplot_tight(2,1,1, [0.04, 0.03]);
    errorbar(v_re_mean(10:end-10,:), v_re_std(10:end-10,:))
    title('Real', 'FontSize', 20);
    xlabel('Time, miliseconds', 'FontSize', 22);
    xticks(ticks);
    xticklabels(label_names);
    legend('comp1','comp2','comp3','comp4');
    pbaspect([5,1,1]);
    set(gca, 'FontSize', 20);
    export_fig([save_dir, figname, '_re.png'], '-png', '-transparent');
    close(gcf);

    errorbar(v_im_mean(10:end-10,:), v_im_std(10:end-10,:))
    title('Imaginary', 'FontSize', 20);
    xlabel('Time, miliseconds', 'FontSize', 22);
    xticks(ticks);
    xticklabels(label_names);
    legend('comp1','comp2','comp3','comp4');
    pbaspect([5,1,1]);
    set(gca, 'FontSize', 20);
    export_fig([save_dir, figname, '_im.png'], '-png', '-transparent');
    close(gcf);
end

% for i = 1:100
%     vv_fl = flip_components(v_re,i);
%     cc = corrcoef(vv_fl);
%     mcc(i) = mean(cc(:));
% end

% cc_ini = corrcoef(v_re);
% [val,key] = max(mcc);


