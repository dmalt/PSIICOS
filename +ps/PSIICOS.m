function [corr, Cpvec, Upwr] = PSIICOS(C, G2dU, SL_rnk,...
                                       sig_rnk, Upwr,...
                                       seed_ind, cp_part,...
                                       is_fast)
% ------------------------------------------------------------------------------------------
% Project from VC and do thresholding on correlations of sources with the cross-spectrum
% ------------------------------------------------------------------------------------------
% FORMAT:
%   [corr, Cpvec, Upwr] = PSIICOS(C, G2dU, SL_rnk, sig_rnk, Upwr, seed_ind, cp_part, is_fast)
% INPUTS:
%   C        - {N_sensors_reduced * N_sensors_reduced x n_Times}
%              sensor-space cross-spectral matrix
%   G2dU     - {N_sensors_reduced x N_sources} forward model matrix 
%              such that each source is served by two columns
%              of this matrix corresponding to the topographies of dipoles
%              in the tangential plane
%   SL_rnk      - scalar; rank of signal leakage subspace. The bigger this value
%              the more data will be removed by the projection from SL. On the
%              contrary, the smaller it is the more SL-related activity will
%              remain in the data. Recommended values are between 350 and 500
%   sig_rnk  - scalar; number of components left after dimensionality reduction
%              of signal subspace
%              if sig_rnk = 0, use mean of cross-spectrum across the time domain 
%   Upwr     - {N_sensors_reduced ^ 2 x SL_rnk} SL subspace basis matrix. 
%              Columns of Upwr span the VC subspace
%   seed_ind - 
%   cp_part  - 
%   is_fast  - 
% OUTPUTS:
%   corr.data  - {1 x Nsrc * (Nsrc - 1) / 2} vector of correlation between
%                topographies and signal subspace
%   corr.IND   -
%   Cpvec      - projected away from the VC subspace sensor space cross-spectral
%                matrix
%   Upwr       - VC subspace basis matrix. Columns of Upwr span the VC subspace
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    %% Preparatory steps
    import ps.ProjectAwayFromPowerComplete

    % --- set defaults --- %
    if(nargin < 8)
        is_fast = true;
    end

    if(nargin < 7)
        cp_part = 'full';
    end

    if(nargin < 6)
        seed_ind = [];
    end

    if(nargin < 5)
        Upwr  = [];
    end;

    if nargin < 4
        sig_rnk = 0;
    end

    if(isempty(SL_rnk))
        if(isempty(Upwr))
            SL_rnk = 350;
        else
            SL_rnk = size(Upwr, 2);
        end;
    end;
    % ------------------- %


    Nsrc = size(G2dU, 2) / 2; % two topography columns per each source of the grid
    n_sensors_C = sqrt(size(C,1));
    n_sensors_G = size(G2dU,1);

    % ---------------- check if matrix sizes are ok. --------------- %
    assert(  n_sensors_C == fix(n_sensors_C),...
            ['NONSQUARE NUMBER OF ROWS IN C: size(C) = ', num2str(size(C))] ); 

    assert(n_sensors_C == n_sensors_G,...
           ['INCONSISTENT NUMBER OF SENSORS IN C AND G2dU: ',...
             num2str(n_sensors_C), ' vs ', num2str(n_sensors_G)]);
    % -------------------------------------------------------------- %

    % -- project coherence matrix from signal-leakage subspace -- %
    if(isempty(Upwr))
        [Cpvec, Upwr] = ProjectAwayFromPowerComplete(C, G2dU, SL_rnk);
    else % use the existing matrix if profided
        assert(n_sensors_C ^ 2 == size(Upwr,1),...
               ['INCONSISTENT SIZES: size(C,1) = ',...
                num2str(size(C,1)), ', size(Upwr,1) = ', num2str(size(Upwr,1))]);
        c = Upwr' * C;
        Cpvec  = C - Upwr * c;
    end;
    % ------------------------------------------------------------ %


    if sig_rnk
        [uc,~,~] = svd(Cpvec, 'econ');
        Cp = uc(:, 1:sig_rnk);
    elseif sig_rnk == 0
        Cp = sum(Cpvec, 2);
        Cp = Cp / norm(Cp, 'fro');
    else
        fprintf('ERROR: PSIICOS: Signal space rank %f\n is not valid', sig_rnk);
        return;
    end

    % --------------------- normalize forward matrix --------------------- %
    for i = 1:Nsrc
        range_i = i * 2 - 1 : i *  2;
        G2dU(:, range_i(1)) = G2dU(:, range_i(1)) / norm(G2dU(:, range_i(1)));
        G2dU(:, range_i(2)) = G2dU(:, range_i(2)) / norm(G2dU(:, range_i(2)));
    end;
    % -------------------------------------------------------------------- %

    if is_fast
        import ps.PSIICOS_ScanFast
        if strcmp(cp_part, 'real')
            [corr.data, corr.IND] = PSIICOS_ScanFast(G2dU, real(Cp));
        elseif strcmp(cp_part, 'imag')
            [corr.data, corr.IND] = PSIICOS_ScanFast(G2dU, imag(Cp));
        elseif strcmp(cp_part, 'full')
            [corr.data, corr.IND] = PSIICOS_ScanFast(G2dU, (Cp));
        else
            fprintf('ERROR: PSIICOS: Option %s for cp_part argument is not valid', sig_rnk);
            return;
        end

        if ~isempty(seed_ind)
            % Add 1 as a placeholder for node coherence with itself
            seed_indices = corr.IND(:,1) == seed_ind | corr.IND(:,2) == seed_ind;
            corr.data = corr.data(seed_indices);
            corr.data = [corr.data(1:seed_ind); 0; corr.data(seed_ind + 1 : end)];
            corr.IND = corr.IND(seed_indices,:);
            corr.IND = [corr.IND(1:seed_ind, :); [seed_ind, seed_ind]; corr.IND(seed_ind + 1 : end, :)];
        end
    else
        [corr.data, corr.IND] = honest_corrs(G2dU, Upwr, Cp, seed_ind, cp_part);
    end

end


function [Cs, IND] = honest_corrs(G2dU, Upwr, Cp, seed_ind, cp_part)
% -------------------------------------- %
% Compute subspace correlations
% -------------------------------------- %

    import ups.indUpperDiag2mat
    n_src = size(G2dU, 2) / 2; % two topography columns per each source of the grid

    if isempty(seed_ind)
        Cs = zeros((n_src ^ 2 - n_src) / 2, 1);
        fprintf('Calculating all-to-all subspace correlations... \n');
        fprintf('Source pair number (max %d): ', (n_src ^ 2 - n_src) / 2);
    else
        Cs = zeros(n_src, 1);
        fprintf('Calculating seed-based subspace correlations... \n');
        fprintf('Source number (max %d): ', n_src);
    end

    p = 1;
    for i_src = 1:n_src
        if isempty(seed_ind)
            for j_src = i_src + 1 : n_src
                Cs(p) = get_ij_subcorr(G2dU, Upwr, Cp, i_src, j_src, cp_part);

                if p > 1
                     for j=0 : log10(p - 1)
                         fprintf('\b');
                     end
                end
                fprintf('%d', p);

                p = p + 1;
            end
        else
            if i_src == seed_ind
                Cs(i_src) = 0;
            else
                Cs(i_src) = get_ij_subcorr(G2dU, Upwr, Cp, i_src, seed_ind, cp_part);
            end

            if i_src > 1
                 for j=0 : log10(i_src - 1)
                     fprintf('\b');
                 end
            end
            fprintf('%d', i_src);
        end
    end
    fprintf(' -> Done\n');

    IND = indUpperDiag2mat(n_src);
    if ~isempty(seed_ind)
        % 
        seed_indices = IND(:,1) == seed_ind | IND(:,2) == seed_ind;
        IND = IND(seed_indices, :);
        IND = [IND(1:seed_ind, :); [seed_ind, seed_ind]; IND(seed_ind + 1 : end, :)];
    end
end


function cs = get_ij_subcorr(G2dU, Upwr, Cp, i_src, j_src, cp_part)
% --------------------------------------------------------------- %
% Get subspace correlation for [i_src, j_src] network
% --------------------------------------------------------------- %
    import ps.subcorr

    assert(i_src ~= j_src,...
          ['ERROR: get_ij_subcorr: i_src = j_src for i_src = ', num2str(i_src),...
          'j_src = ', num2str(j_src)])

    ai = G2dU(:, i_src * 2 - 1 : i_src * 2);
    aj = G2dU(:, j_src * 2 - 1 : j_src * 2);
    Gij = kron(ai, aj);
    Gji = kron(aj, ai);

    if strcmp(cp_part,'real')
        G_re = Gij + Gji;
        temp = Upwr' * G_re;
        G_re = G_re - Upwr * temp;
        subc = subcorr(G_re, real(Cp));
    elseif strcmp(cp_part, 'imag')
        G_im = Gij - Gji;
        subc = subcorr(G_im, imag(Cp));
    elseif strcmp(cp_part, 'full')
        G_re = Gij + Gji;
        G_re = G_re - Upwr * (Upwr' * G_re);
        subc_re = subcorr(G_re, real(Cp));
        G_im = Gij - Gji;
        subc_im = subcorr(G_im, imag(Cp));
        subc = sqrt(subc_re .^ 2 + subc_im .^ 2);
    else
        error(['Unknown option ', cp_part]);
    end

    cs = subc(1);
end
