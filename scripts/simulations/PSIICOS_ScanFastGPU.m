function [Cs, IND] = PSIICOS_ScanFastGPU(G2dU, Cp, is_imag, blocksize)
% -------------------------------------------------------------------------
% Perform scanning algorithm to find strongest connections using correlation
% of cross-spectrum with the forward operator
% -------------------------------------------------------------------------
% FORMAT:
%   [Cs, IND Cs0] = PSIICOS_ScanFast(G2dU, Cp)
% INPUTS:
%   G2dU       - {n_sensors x n_sources} matrix of forward model
%   Cp         - {n_sensors ^ 2 x n_times} or
%                {n_sensors ^ 2 x n_components}
%                matrix of timeseries or left singular vectors of
%                cross-spectrum on sensors
% OUTPUTS:
%   Cs         - {(n_sources ^ 2 - n_sources) / 2 x 1} vector of
%                correlations between source topographies
%                and forward operator
%   IND        - {(n_sources ^ 2 - n_sources) / 2 x 2} matrix of
%                indices to build a mapping between upper
%                triangle and (i,j) matrix indexing
%                IND(l,:) --> [i,j]
% ___________________________________________________________________________
% Alex Ossadtchii, ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    if nargin < 3
        is_imag = false;
    end

    [Nsns, Nsrc2] = size(G2dU);
    Nsrc = Nsrc2 / 2;

    % if(size(Cp,1)~= Nsns)
    %     disp('Incomptible dimensions G2dU vs Cp');
    %     return;
    % end

    n_comp = size(Cp, 2);
    T = zeros(n_comp, Nsrc * (Nsrc - 1) / 2);
    D = zeros(n_comp, Nsrc * (Nsrc - 1) / 2);
    IND = zeros(Nsrc * (Nsrc - 1) / 2, 2);

    % below is the optimized implementation of this:
    % Look at each pair and estimate subspace correlation
    % between cross-spectrum and topography of this pair
    % tic
    % p = 1;
    % for i=1:Nsrc
    %     range_i = i * 2 - 1 : i * 2;
    %     ai = G2dU(:,range_i)';
    %     for j=i + 1:Nsrc
    %          range_j = j * 2 - 1 : j * 2;
    %          aj = G2dU(:, range_j)';
    %         cs = ai * Cp * aj';
    %         [u s v] = svd(cs);
    %         Cs0(p) = max(diag(s));
    %          p = p + 1;
    %     end;
    % end;
    % toc
    %
    G2dU_gpu = gpuArray(single(G2dU));
    % G2dU_gpu = (single(G2dU));
    Cs = zeros(Nsrc, 'single', 'gpuArray');
    blocksize_x2 = blocksize * 2;
    n_blocks = floor(Nsrc / blocksize); % number of equal blocks; remaining elements are treated separately
    is_last_block = mod(Nsrc, blocksize);



    % p = 1;
    for iComp = 1:n_comp
        for ii = 1:n_blocks
            Cp_sq = gpuArray(single(reshape(Cp(:,iComp), Nsns, Nsns)));
            % Cp_sq =(single(reshape(Cp(:,iComp), Nsns, Nsns)));
            % ss = G2dU_gpu(:,1 + blocksize * (ii - 1): blocksize * ii)' * Cp_sq * G2dU_gpu;

            a11 = G2dU_gpu(:,1 + blocksize_x2 * (ii - 1):2: blocksize_x2 * ii)' * Cp_sq * G2dU_gpu(:,1:2:end);
            a12 = G2dU_gpu(:,1 + blocksize_x2 * (ii - 1):2: blocksize_x2 * ii)' * Cp_sq * G2dU_gpu(:,2:2:end);
            a21 = G2dU_gpu(:,2 + blocksize_x2 * (ii - 1):2: blocksize_x2 * ii)' * Cp_sq * G2dU_gpu(:,1:2:end);
            a22 = G2dU_gpu(:,2 + blocksize_x2 * (ii - 1):2: blocksize_x2 * ii)' * Cp_sq * G2dU_gpu(:,2:2:end);

            a11_c = conj(a11);
            a22_c = conj(a22);
            a12_c = conj(a12);
            a21_c = conj(a21);
            % clear ss;

            a11_2 = real(a11 .* a11_c);
            a22_2 = real(a22 .* a22_c);
            a12_2 = real(a12 .* a12_c);
            a21_2 = real(a21 .* a21_c);

            ll = 2 * real(a11 .* a22 .* a12_c .* a21_c);

            clear a11;
            clear a12;
            clear a22;
            clear a21;
            clear a11_c;
            clear a22_c;
            clear a12_c;
            clear a21_c;

            % Cs = abs(0.5 * (a11 + a22) + sqrt(0.25 * (a11 - a22) .^ 2 + a12 .* a21)) .^ 2;
            Cs(1 + blocksize * (ii - 1): blocksize * ii, :) = (0.5 * (a11_2 + a12_2 + a21_2 + a22_2) + ...
            sqrt( 0.25 * (a11_2 + a12_2 - a21_2 - a22_2) .^ 2 +  a11_2 .* a21_2 + ll + a22_2 .* a12_2));

            clear a11;
            clear a22;
            clear a12;
            clear a21;
            % ss = sparce(triu(ss,1));
            % ss_c = conj(ss);
            % ss_long = real(ss .* ss_c);
            % ss_longd = ss(1:2:end,:) .* (ss_c(2:2:end,:));
            % clear ss;
            % c2_11_22 = ss_long(:,1:2:end) + ss_long(:,2:2:end);
            % c2_12_21 = ss_longd(:,1:2:end) + ss_longd(:,2:2:end);
            % clear ss_long;
            % clear ss_longd;
            % c2_12_21 = real(c2_12_21 .* conj(c2_12_21));
            % T = c2_11_22(1:2:end,:) + c2_11_22(2:2:end,:);
            % D = c2_11_22(1:2:end,:) .* c2_11_22(2:2:end,:) - c2_12_21;
            % clear c2_11_22;
            % clear c2_12_21;
            % Cs(1 + blocksize_2 * (ii - 1): blocksize_2 * ii, :) = sqrt(0.5 * T + sqrt(0.25 * T .* T - D));
            % % Cs =
            % clear T
            % clear D
        end
        % -------------------- Deal with the remaining elements -------------------- %
        if is_last_block
            Cp_sq = gpuArray(single(reshape(Cp(:,iComp), Nsns, Nsns)));
            % Cp_sq =(single(reshape(Cp(:,iComp), Nsns, Nsns)));
            % ss = G2dU_gpu(:,1 + blocksize * (ii - 1): blocksize * ii)' * Cp_sq * G2dU_gpu;

            a11 = G2dU_gpu(:,1 + blocksize_x2 * n_blocks:2: end)' * Cp_sq * G2dU_gpu(:,1:2:end);
            a12 = G2dU_gpu(:,1 + blocksize_x2 * n_blocks:2: end)' * Cp_sq * G2dU_gpu(:,2:2:end);
            a21 = G2dU_gpu(:,2 + blocksize_x2 * n_blocks:2: end)' * Cp_sq * G2dU_gpu(:,1:2:end);
            a22 = G2dU_gpu(:,2 + blocksize_x2 * n_blocks:2: end)' * Cp_sq * G2dU_gpu(:,2:2:end);

            a11_c = conj(a11);
            a22_c = conj(a22);
            a12_c = conj(a12);
            a21_c = conj(a21);
            % clear ss;

            a11_2 = real(a11 .* a11_c);
            a22_2 = real(a22 .* a22_c);
            a12_2 = real(a12 .* a12_c);
            a21_2 = real(a21 .* a21_c);

            ll = 2 * real(a11 .* a22 .* a12_c .* a21_c);

            clear a11;
            clear a12;
            clear a22;
            clear a21;
            clear a11_c;
            clear a22_c;
            clear a12_c;
            clear a21_c;

            % Cs = abs(0.5 * (a11 + a22) + sqrt(0.25 * (a11 - a22) .^ 2 + a12 .* a21)) .^ 2;
            Cs(1 + blocksize * n_blocks:end, :) = (0.5 * (a11_2 + a12_2 + a21_2 + a22_2) + ...
            sqrt( 0.25 * (a11_2 + a12_2 - a21_2 - a22_2) .^ 2 +  a11_2 .* a21_2 + ll + a22_2 .* a12_2));

            clear a11;
            clear a22;
            clear a12;
            clear a21;
        end
        % -------------------------------------------------------------------------- %

        Cs = gather(Cs);
        Cs = nonzeros(triu(Cs,1)');


        % for iSrc = 1:Nsrc
        %     % ---- Take iSrc-th location topographies ---- %
        %     % ai = G2dU(:, iSrc * 2 - 1 : iSrc * 2)';
        %     % tmp = ai * Cp_sq;
        %     % cslong = tmp * G2dU;
        %     % cslong = ss(iSrc*2-1:iSrc*2,:);
        %     % if is_imag
        %     %     cslong = imag(cslong);
        %     % end
        %     % cs2long = cslong .* conj(cslong);
        %     cs2long = ss_long(iSrc*2-1:iSrc*2,:);
        %     % cs2longd = cslong(1,:) .* conj(cslong(2,:));
        %     cs2longd = ss_longd(iSrc,:);
        %     cs2_11_22 = [sum(reshape(cs2long(1,:), 2, Nsrc), 1);...
        %                  sum(reshape(cs2long(2,:), 2, Nsrc), 1)];

        %     cs2_12_21 = sum(reshape(cs2longd, 2, Nsrc), 1);
        %     Ti = sum(cs2_11_22, 1);
        %     Di = prod(cs2_11_22, 1) - cs2_12_21 .* conj(cs2_12_21);
        %     T(iComp, p : p + Nsrc - 1 - iSrc) = Ti(iSrc + 1 : Nsrc);
        %     D(iComp, p : p + Nsrc - 1 - iSrc) = Di(iSrc + 1 : Nsrc);
        % end
        % IND(p : p + Nsrc - iSrc - 1, 2) = (iSrc + 1 : Nsrc)';
        % IND(p : p + Nsrc - iSrc - 1, 1) = iSrc;
        % p = p + Nsrc - iSrc;
    end
    % Cs = 0.5 * T + sqrt(0.25 * T .* T - D);
    % Cs = sum(Cs, 1);
    % Cs = max(Cs, [], 1)';
end
