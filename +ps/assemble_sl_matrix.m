function S = assemble_sl_matrix(G, normalize, loose)
    n_src = size(G, 2); 
    n_sen = size(G, 1);
    % project for each potential source
    fprintf('Collecting signal leakage subspace...\n');
    if loose
        n_src = n_src / 2; % in loose case we have two topography columns per source
        S = zeros(n_sen ^ 2, n_src * 3);
        for i = 1:n_src
             gi = G(:, 2 * i - 1);
             gj = G(:, 2 * i);

             S(:, 3 * i - 2) = kron(gi, gi);
             S(:, 3 * i - 1) = kron(gj, gj);
             S(:, 3 * i) = kron(gi, gj) + kron(gj, gi);
        end
    else
        S = zeros(n_sen ^ 2, n_src);
        for i = 1 : n_src
             gi = G(:, i);
             S(:, i) = kron(gi, gi);
        end
    end

    if normalize
        S = S ./ sum(S .* 2, 1); % normalize each column
    end
end
