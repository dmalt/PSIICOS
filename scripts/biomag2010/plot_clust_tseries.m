function plot_clust_tseries(con_clust, con_all, HM, CTs, band_name, part, i_comp)
% ----------------------------------------------------------------------- %
% Plot errorbar timeseries for first 3 clusters in clusterized connections
% ----------------------------------------------------------------------- %

    % -------------- define ticks ------------- %
    sfreq = 1200;
    ticks = 0:150:1200;
    label_vals = ticks / sfreq * 1000;
    label_names = cell(size(label_vals));
    for ii = 1:length(label_vals)
        label_names{ii} = num2str(label_vals(ii));
    end
    % ----------------------------------------- %
    colorscheme = ups.plt.GetColors(3);
    colors = {colorscheme(2:end).RGB};
    all_con = con_all.Merge().conInds{1};
    G = HM.gain;
    n_sen = size(G,1);

    n_clust = [];



    for i_clust = 1:length(con_clust.conInds)
        if length(con_clust.conInds{i_clust}) > 15
            n_clust = [n_clust,i_clust];
        end
    end
    figname = ['dics_tseries_', band_name, '_', part, '_', num2str(i_comp)]
    if ~isempty(n_clust)
        figure('Name', figname);
    end

    for kk = 1:length(n_clust)
    % for i_clust = 1:3
        i_clust = n_clust(kk);
        clust_inds = con_clust.conInds{i_clust};
        i_clust_col = con_clust(i_clust).icol;

        for i_stick = 1:size(clust_inds, 1)
            lr = sum(clust_inds(i_stick,:) == all_con, 2) == 2;
            rl = sum(clust_inds(i_stick,:) == all_con(:,[2,1]), 2) == 2;
            inds = find(lr | rl);
            for jj = 1:length(inds)
                ind(i_stick) = inds(jj);
            end
        end


        n_clust_con = size(clust_inds,1);
        ts = zeros(n_clust_con, size(CTs{1},2));

        for i_resamp = 1:n_clust_con
            c = all_con(ind(i_resamp),:);

            ii = c(1);
            jj = c(2);

            Gi = G(:, 2 * ii - 1:2 * ii);
            Gj = G(:, 2 * jj - 1:2 * jj);

            CTk = mean(CTs{ind(i_resamp)},2);

            CC = reshape(CTk, n_sen, n_sen);

            [A,Or] = ups.conn.DICS(CC, [Gi,Gj]);

            Ai = Or{1}' * A(1:2,:);
            Aj = Or{2}' * A(3:4,:);
            AA = kron(Ai,Aj);

            ts(i_resamp,:) = AA * CTs{ind(i_resamp)};
        end

        ats = abs(ts);
        mts = mean(ats,1);
        stdts = std(ats,[],1) / sqrt(n_clust_con);

        icol_ = mod(i_clust_col, length(colors));

        if icol_ == 0
           icol_ = length(colors);
        end

        hold on;
        h = errorbar(mts(10:end-10)',stdts(10:end-10)');
        set(h, 'Color',  colors{icol_})%, 'LineWidth', 1);
    end

    if ~isempty(n_clust)
        title(part, 'FontSize', 20);
        xticks(ticks);
        xticklabels(label_names);
        % legend('clust1','clust2','clust3');
        xlabel('Time, miliseconds', 'FontSize', 22);
        pbaspect([5,1,1]);
        set(gca, 'FontSize', 20);
        export_fig([figname, '.png'], '-png', '-transparent');
        close(gcf);
    end
end
