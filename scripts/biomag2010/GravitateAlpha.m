function av_bundles = GravitateAlpha(bundles, c_dist, mass)
% ------------------------------------------------------------ %
% For each stick in bundles get average and set alpha based on
% proximity to other sticks
% ------------------------------------------------------------ %
    % av_bundles = bundles;
    base_alpha = 0.2;
    size_inc_alpha = 0.02;
    G = 3e-8;
    n_bundles = length(bundles.conInds);
    mass = zeros(n_bundles,1);
    stick_ends = cell(n_bundles,1);
    for i_b = 1:n_bundles
        n_sticks = size(bundles(i_b).conInds{1},1);
        mass(i_b) = n_sticks;
        stick = bundles(i_b).Average();
        stick_ends{i_b} = [stick.headModel.GridLoc(stick.conInds{1}(1),:);...
                           stick.headModel.GridLoc(stick.conInds{1}(2),:)];
    end
    av_bundles = bundles.Average();
    n_ends = length(stick_ends);
    dist_matr = zeros(n_ends);
    for ii = 1:n_ends
        for jj = 1:n_ends
            dist_matr(ii,jj) = StickDist(stick_ends{ii}, stick_ends{jj});
        end
    end
    dist_matr = dist_matr + eye(size(dist_matr));
    AA = ones(size(dist_matr)) ./ (dist_matr .^ 2);
    AA = AA - diag(AA);

    influence = zeros(n_ends,1);
    influence = AA * mass * G + base_alpha + eye(size(AA)) * mass * size_inc_alpha;
    % ind_min = influence < eye(size(AA)) * mass * base_alpha;
    % influence(ind_min)= mass(ind_min) * base_alpha;
    influence(influence > 1) = 1;

    av_bundles.alpha = influence';


end
