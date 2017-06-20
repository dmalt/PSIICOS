c_mask_inc = data_hr > clust_threshold;
c_mask_dec = data_hr < -clust_threshold;

[c_inc, c_dec] = get_clusters(mat, c_mask_inc, c_mask_dec);


stat_inc_true = zeros(length(c_inc), 1);
for iclust = 1:length(c_inc)
    stat_inc_true(iclust) = sum(data_hr(c_inc{iclust}));
end

stat_dec_true = zeros(length(c_dec), 1);
for iclust = 1:length(c_dec)
    stat_dec_true(iclust) = abs(sum(data_hr(c_dec{iclust})));
end
max_stat_true = sort([stat_inc_true; stat_dec_true], 'descend');

% mask = ones(2004, 1);
% mask = ones(15002, 1);
% cl_num = 5;
% mask(c_inc{cl_num}) = zeros(length(c_inc{cl_num}), 1);

% figure;
% h = plot_brain_cmap_hemisplit(CtxHHR_dst, CtxHR_dst, [], data_hr,...
%                               mask, 0.1, seed_dst_xyz_approx);

figure;
h = plot_brain_cmap_hemisplit(CtxHHR_dst, CtxHR_dst, [], data_hr,...
                              ~c_mask_inc, 0.1, seed_dst_xyz_approx);
