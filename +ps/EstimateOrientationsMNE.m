function [G, Gw, theta, lambda_ratio, pwr] = EstimateOrientationsMNE(G2, C, lambda)
% Estimate oriented version of G2 using minimum norm filters
%
    D = trace(G2 * G2') / size(G2, 1);
    Gamma_reg = G2 * G2' + lambda * D * eye(size(G2, 1));
    % fprintf('condition number: %.2f\n', cond(Gamma_reg))
    C_src_x_G2 = (Gamma_reg \ C / Gamma_reg) * G2;
    n_sen = size(G2, 1);
    n_src = size(G2, 2) / 2;
    G = zeros(n_sen, n_src);
    Gw = zeros(n_sen, n_src);
    if nargout > 2
        theta = zeros(n_src, 1);
    end
    if nargout > 3
        lambda_ratio = zeros(n_src, 1);
    end
    if nargout > 4
        pwr = zeros(n_src, 1);
    end
    for i = 1:n_src
        gi = G2(:, 2 * i - 1 : 2 * i);
        ti = C_src_x_G2(:, 2 * i - 1 : 2 * i);
        m = gi' * ti;
        [v, d] = eig(m);
        [~, ind] = sort(diag(d), 'descend');
        ori = v(:, ind(1));
        G(:, i) = gi * ori;% * d(1);
        Gw(:, i) = gi * ori * sqrt(d(ind(1), ind(1)));
        if nargout > 2
            theta(i) = v(1, ind(1)) / v(2, ind(1));
        end
        if nargout > 3
            lambda_ratio(i) = sqrt(d(ind(2), ind(2)) / d(ind(1), ind(1)));
        end
        if nargout > 4
            pwr(i) = d(ind(1), ind(1));
        end
    end
end
