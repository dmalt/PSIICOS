function corr = PSIICOS_normalized(Trials, G2, lambda)
    n_tr = size(Trials, 3);
    n_times = size(Trials, 2);
    is_induced = false;
    n_src = size(G2, 2) / 2;
    n_sen = size(G2, 1);
    count = 0;

    IND = ups.indUpperDiag2mat(n_src);
    % PLV = zeros(size(IND, 1), 1);
    % pwr = gpuArray(zeros(size(IND, 1), 1));
    % cp = gpuArray(zeros(size(IND, 1), 1));
    pwr = gpuArray(zeros(n_src, n_src, 'single'));
    cp_re = gpuArray(zeros(n_src, n_src, 'single'));
    cp_im = gpuArray(zeros(n_src, n_src, 'single'));
    CT_temp = ups.conn.CrossSpectralTimeseries(Trials, is_induced);
    CC = reshape(mean(CT_temp, 2), n_sen, n_sen);
    [corr, U, S, G, W_re, W_im] = ps.PSIICOS_LORETA(CC, G2, lambda);
    tmp = gpuArray(single((U * S) * U'));
    index = tril(true([n_src, n_src]), -1);
    W_re = gpuArray(single(W_re));
    W_im = gpuArray(single(W_im));
    G = gpuArray(single(G));
    for i = 1:n_tr - 1
        % fprintf('Trial - > %d\n', i);
        CT = gpuArray(single(...
            ups.conn.CrossSpectralTimeseries(Trials(:, :, i:i+1),...
                is_induced, false)));
        % a = tic();
        % fprintf('Projection\n');
        CT_proj = reshape(tmp * CT, n_sen, n_sen, n_times);
        % fprintf('Elapsed %.2f\n', toc(a))
        b = tic;
        for t = 1:n_times
            % C = reshape(CT_proj(:, t), n_sen, n_sen);
            C = CT_proj(:, :, t);
            % if i == 1 && t == 1
            %     [corr, U, S, G] = ps.PSIICOS_LORETA(C, G2, lambda);
            % else
            % [corr] = ps.PSIICOS_LORETA(C, G2, lambda, U, S, G, W_re, W_im);
            % C_nosl = reshape(C_vec_nosl, n_sen, n_sen);
            % C_nosl = imag(C);
            % aa = tic;
            C_src_re = G' * (real(C) * G);
            % C_src_re = C_src_re(index);
            C_src_im = G' * (imag(C) * G);
            % C_src_im = C_src_im(index);
            % fprintf('aa %.2f\n', toc(aa) * 1000)
            % corr.data = C_src(i) ./ W(i);
            % corr.data = real(C_src(i)) ./ W_re(i) + 1i * imag(C_src(i)) ./ (W_im(i).^2) .* W_re(i);
            % corr.data = C_src(i) ./ sqrt(W_re(i).^2 + W_im(i) .^ 2);
            % corr.data = C_src(i) ./ sqrt(W_re(i).^2 + W_im(i) .^ 2);
            % bb = tic;
            c_re = C_src_re ./ W_re;
            c_im = C_src_im ./ W_im;
            % fprintf('bb %.2f\n', toc(bb) * 1000)
            % end
            % c = corr.data;
            % disp(c(1:20))
            % cc = tic;
            %
            % PLV = PLV + c ./ abs_c;
            pwr = pwr + sqrt(c_re .^ 2 + c_im .^ 2);
            cp_re = cp_re + c_re;
            cp_im = cp_im + c_im;
            % fprintf('cc %.2f\n', toc(cc) * 1000)
            % fprintf('%d-->%d\n', i, t);
            count = count + 1;
        end
        % fprintf('Inner loop %.2f\n', toc(b))
    end
    % PLV = PLV / count;
    % corr.PLV = PLV;
    cp = cp_re + 1i* cp_im;
    corr.cp = gather(cp(index));
    corr.pwr = gather(pwr(index));
    corr.data = corr.cp ./ corr.pwr;
    corr.IND = IND;
end

function C_src = compute_source_C(C, G, SS_reg_inv, U, n_sen)
    C_vec = (C(:));
    C_vec_nosl = U * SS_reg_inv * (U' * C_vec);
    C_nosl = reshape(C_vec_nosl, n_sen, n_sen);
    % C_nosl = imag(C);
    C_src = G' * (C_nosl * G);
end

