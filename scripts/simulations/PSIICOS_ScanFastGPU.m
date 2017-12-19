function [Cs, IND] = PSIICOS_ScanFastGPU(G2dU, Cp, is_imag)
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
    Cs = zeros(1, Nsrc * (Nsrc - 1) / 2, 'single', 'gpuArray');
    blocksize = 2000;
    blocksize_2 = blocksize / 2;


    % p = 1;
    jj = 1;
    for iComp = 1:n_comp
        for ii = 1:15
            Cp_sq = gpuArray(single(reshape(Cp(:,iComp), Nsns, Nsns)));
            % Cp_sq =(single(reshape(Cp(:,iComp), Nsns, Nsns)));
            ss = G2dU_gpu(:, 1 + blocksize * (ii - 1): blocksize * ii)' * Cp_sq * G2dU_gpu;
            ss_long = ss .* conj(ss);
            ss_longd = ss(1:2:end,:) .* conj(ss(2:2:end,:));
            clear ss;
            c2_11_22 = ss_long(:,1:2:end) + ss_long(:,2:2:end);
            c2_12_21 = ss_longd(:,1:2:end) + ss_longd(:,2:2:end);
            clear ss_long;
            clear ss_longd;
            c2_12_21 = c2_12_21 .* conj(c2_12_21);

            T = c2_11_22(1:2:end,:) + c2_11_22(2:2:end,:);
            mask = triu(true(size(T)), 1 + blocksize_2 * (ii - 1))';
            T = T';
            T = T(mask);
            ll = length(T);

            D = c2_11_22(1:2:end,:) .* c2_11_22(2:2:end,:) - c2_12_21;
            D = D';
            D = D(mask);

            clear c2_11_22;
            clear c2_12_21;
            Cs(jj:jj + ll - 1) = 0.5 * T + sqrt(0.25 * T .* T - D);
            jj = jj + ll;
            % Cs =
            clear T
            clear D
        end

        Cs = gather(Cs)';
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
