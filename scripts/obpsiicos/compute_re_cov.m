function C_sl = compute_re_cov(G)

    rank_g = rank(G);
    C_sl_gpu =  zeros(rank_g ^ 2, 'gpuArray');
    % C_sl =  zeros(rank_g ^ 2);
    n_sites = length(G) / 2;
    % n_ch = size(G, 1);
    Swp_opr = zeros(4);
    Swp_opr(1,1) = 1;
    Swp_opr(4,4) = 1;
    Swp_opr(2,3) = 1;
    Swp_opr(3,2) = 1;

    for i = 1:n_sites
         range_i = i * 2 - 1 : i * 2;
         ai = single(G(:, range_i));
         rng = 1:4;
         % XX_gpu = zeros(rank_g ^ 2, 4 * (n_sites - i), 'single', 'gpuArray');
         XX = zeros(rank_g ^ 2, 4 * (n_sites - i), 'single');
         for j = i + 1 : n_sites
            range_j = j * 2 - 1 : j * 2;
            aj = single(G(:, range_j));
            XX(:, rng) = kron(ai, aj) + kron(aj, ai) * Swp_opr;
            % rank(XX)
            rng = rng + 4;
        end
        % rank(XX)
        % size(XX)
        XX_gpu = gpuArray(XX);
        C_sl_gpu = C_sl_gpu + (XX_gpu * XX_gpu');
        % C_sl = C_sl + (XX * XX');
        % rank(C_sl)
        disp(i)
    end
    C_sl = gather(C_sl_gpu);
end
