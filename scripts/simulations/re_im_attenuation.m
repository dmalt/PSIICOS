% ---------------------------------------------------------------- %
% See how real projected and imaginary 2-topographies decrease with
% the increase of correlation between topos of each node in a pair
% ---------------------------------------------------------------- %
% Date: 2017-12-14
% Author: dmalt
% ________________________________________________________________ %

mem_projector =  memoize(@ps.GetSLProjectorMegLoose);
mem_dim_reducer = memoize(@ups.ReduceSensorSpace);
mem_tangent = memoize(@ups.ReduceToTangentSpace);

%  setup paths and load data {{{setup %
colors = ups.plt.GetColors(3);
colors_rgb = {colors(1:end).RGB};

sim_data_path = '~/Code/matlab/utils_psiicos/InputData4Simulations.mat';
ch_path =  '~/Code/matlab/utils_psiicos/channel_vectorview306.mat';

sim_data = load(sim_data_path);
ch = load(ch_path);
chans = ch.Channel;

GainSVDTh = 0.01;
% ch_used = ups.bst.PickChannels(chans, {'MEG GRAD', 'MEG MAG'});
ch_used = ups.bst.PickChannels(chans, 'MEG GRAD');
G = sim_data.G.Gain(ch_used,:);
Gr = mem_dim_reducer(G, GainSVDTh);
G2 = mem_tangent(Gr);

pwr_rnks = {10, 50, 90, 150, 250, 350, 500, 700};

%  setup}}} %
n_monte_carlo = 10000;
n_re = cell(length(pwr_rnks), 1);
cc_re = cell(length(pwr_rnks), 1);
cc_im = zeros(n_monte_carlo, 1);

zeros(n_monte_carlo, 1);
n_im = zeros(n_monte_carlo, 1);
[n_sen, n_dip] = size(G2);

for i_monte = 1:n_monte_carlo
    ij = randi(n_dip, [2,1]);
    gi = G2(:,ij(1));
    gi = gi / norm(gi);
    gj = G2(:,ij(2));
    gj = gj / norm(gj);

    G_ij_im = kron(gi,gj) - kron(gj,gi);
    n_im(i_monte) = norm(G_ij_im);

    % cc_im(i_monte) = abs(corr(gi,gj));
    cc_im(i_monte) = (corr(gi,gj));
    % i_monte
end

for i_rnk = 1:length(pwr_rnks)
    Upwr = ps.GetSLProjectorMegLoose(G2, pwr_rnks{i_rnk});
    n_re{i_rnk} = zeros(n_monte_carlo, 1);
    cc_re{i_rnk} = zeros(n_monte_carlo, 1);

    %  monte-carlo iterations {{{monte %
    for i_monte = 1:n_monte_carlo
        ij = randi(n_dip, [2,1]);
        gi = G2(:,ij(1));
        gi = gi / norm(gi);
        gj = G2(:,ij(2));
        gj = gj / norm(gj);
        G_ij_re = kron(gi, gj) + kron(gj, gi);
        G_ij_re = G_ij_re - Upwr * (Upwr' * G_ij_re);
        n_re{i_rnk}(i_monte) = norm(G_ij_re);
        % cc_re{i_rnk}(i_monte) = abs(corr(gi,gj));
        cc_re{i_rnk}(i_monte) = (corr(gi,gj));
        % i_monte
    end
    %  monte}}} %
end

figure;
for i_rnk = 1:length(pwr_rnks)
    [~,key] = sort(cc_re{i_rnk}, 'descend');
    % figure;
    % scatter(cc(key), n_re(key), 'filled', 'MarkerFaceColor', colors_rgb{2}, 'MarkerEdgeColor', colors_rgb{1});
    % hold on;
    % scatter(cc(key), n_im(key), 'filled', 'MarkerFaceColor', colors_rgb{3}, 'MarkerEdgeColor', colors_rgb{1});
    p = polyfit(cc_re{i_rnk}(key),n_re{i_rnk}(key), 2);
    pp = polyval(p, cc_re{i_rnk}(key));
    hold on;
    h = plot(cc_re{i_rnk}(key), pp);
    % i_
    % set(h, 'Color', colors_rgb{2}, 'LineWidth', 3);
    set(h, 'LineWidth', 2);
end

[~,key] = sort(cc_im, 'descend');
p = polyfit(cc_im(key),n_im(key), 2);
pp = polyval(p, cc_im(key));
hold on;
h = plot(cc_im(key), pp, '-.');
set(h, 'LineWidth', 2);

pwr_rnk_names = cellfun(@num2str, pwr_rnks, 'UniformOutput', false);
legend(pwr_rnk_names{:}, 'im');

[~,key] = sort(cc_re{4});
hold on;

figure;
scatter(cc_re{6}(key), n_re{6}(key), 'filled', 'MarkerFaceColor', colors_rgb{2}, 'MarkerEdgeColor', colors_rgb{1});
hold on;
scatter(cc_im, n_im, 'filled', 'MarkerFaceColor', colors_rgb{3}, 'MarkerEdgeColor', colors_rgb{1});
