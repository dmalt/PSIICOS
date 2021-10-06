function [corr, U, SS_reg_inv, G, W_re, W_im] = PSIICOS_LORETA(C, G2, lambda, U, SS_reg_inv, G, W_re, W_im)
    C_re = real(C);
    if nargin < 6
        [G, Gw] = ps.EstimateOrientationsMNE(G2, C_re, 1);
    end

    normalize = false;
    loose = false;
    if nargin < 4
        [U, S] = ps.ComputeSlSvd(G, normalize, loose);
        D = trace(S * S') / size(S, 1);
        SS_reg_inv = diag(1 ./ diag(S * S' + lambda * D * eye(size(S, 1))));
        fprintf('condition number: %.2f\n', cond(SS_reg_inv))
    end

    [n_sen, n_src] = size(G);
    if nargin < 7
        [W_re, W_im] = compute_weights(G, SS_reg_inv, U, n_sen, n_src);
    end
    inner = tic;
    C_src = compute_source_C(C, G, SS_reg_inv, U, n_sen);

    i = tril(true(size(C_src)), -1);
    % corr.data = C_src(i) ./ W(i);
    % corr.data = real(C_src(i)) ./ W_re(i) + 1i * imag(C_src(i)) ./ (W_im(i).^2) .* W_re(i);
    % corr.data = C_src(i) ./ sqrt(W_re(i).^2 + W_im(i) .^ 2);
    % corr.data = C_src(i) ./ sqrt(W_re(i).^2 + W_im(i) .^ 2);
    %
    corr.data = real(C_src(i)) ./ W_re(i) + 1i * imag(C_src(i)) ./ W_im(i);
    % corr.data = 1i * imag(C_src(i)) ./ W_im(i);
    corr.IND = ups.indUpperDiag2mat(n_src);
    toc(inner)
end


function [W_re, W_im] = compute_weights(G, SS_reg_inv, U, n_sen, n_src)
    % fprintf('')
    % disp(size(U))
    % disp(size(SS_reg_inv))
    tmp = U * sqrt(SS_reg_inv);
    W_re = zeros(n_src, n_src);
    symm_count = n_sen * (n_sen + 1) / 2;

    for i = 1:symm_count
        Ur = reshape(tmp(:, i), n_sen, n_sen);
        W_re = W_re + (G' * Ur * G) .^ 2;
        % disp(norm(Ur + Ur', 'fro'))
    end
    W_re = sqrt(W_re);

    W_im = zeros(n_src, n_src);
    for i = symm_count + 1 : n_sen ^ 2
        Ur = reshape(tmp(:, i), n_sen, n_sen);


        % disp(norm(Ur + Ur', 'fro'))
        W_im = W_im + (G' * Ur * G) .^ 2;
    end
    W_im = sqrt(W_im);
end


function C_src = compute_source_C(C, G, SS_reg_inv, U, n_sen)
    C_vec = (C(:));
    C_vec_nosl = U * SS_reg_inv * (U' * C_vec);
    C_nosl = reshape(C_vec_nosl, n_sen, n_sen);
    % C_nosl = imag(C);
    C_src = G' * (C_nosl * G);
end

