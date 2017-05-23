function [clusters_inc, clusters_dec] = get_clusters(adj_mat, mask_inc, mask_dec)
% ------------------------------------------------
% Find clusters
% ------------------------------------------------

    [n_src, n_src] = size(adj_mat);
    adj_mat_inc = adj_mat;
    adj_mat_dec = adj_mat;

    adj_mat_inc(~mask_inc,:) = zeros(sum(~mask_inc), n_src);
    adj_mat_inc(:,~mask_inc) = zeros(n_src, sum(~mask_inc));

    adj_mat_dec(~mask_dec,:) = zeros(sum(~mask_dec), n_src);
    adj_mat_dec(:,~mask_dec) = zeros(n_src, sum(~mask_dec));


    clusters_inc = extract_nontrivial_clusters(adj_mat_inc);
    clusters_dec = extract_nontrivial_clusters(adj_mat_dec);
end


function clusters = extract_nontrivial_clusters(adj_mat)
    import ups.conncomp

    clusters = {};
    clust_num = 0;
    [n_clust, C] = conncomp(adj_mat);
    for i_clust = 1:n_clust
        clust_idx = find(C == i_clust);
        if length(cust_idx) > 1
            clust_num = clust_num + 1;
            clusters{clust_num} = clust_idx;
        end
    end
end
