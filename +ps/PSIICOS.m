function [corr, Cpvec, Upwr] = PSIICOS(CT, G2, sl_rnk, sig_rnk, Upwr,...
                                       seed_ind, cp_part, is_fast)
% ----------------------------------------------------------------------
% Project from signal leakage and compute subcors of source-level pairs
% with the cross-spectrum
% ----------------------------------------------------------------------
% FORMAT:
%   [corr, Cpvec, Upwr] = PSIICOS(CT, G2, sl_rnk, sig_rnk, Upwr,...
%                                 seed_ind, cp_part, is_fast)
% INPUTS:
%   CT     - {n_sensors * n_sensors x n_times}
%            sensor-space cross-spectral matrix
%   G2     - {n_sensors x n_sources * 2} matrix;
%            MEG loose orientation forward operator
%            projected to the tangential plane
%   sl_rnk - scalar; rank of signal leakage subspace.
%            The bigger this value the more data will be
%            removed by the projection from SL-subspace but
%            'good' signal in real part will be affected
%            proportionally.
%            Recommended values are between 150 and 500.
%            default=350;
%   sig_rnk  - scalar; number of components left after
%              dimensionality reduction of signal subspace;
%              if sig_rnk = 0, use mean of cross-spectrum
%              across time
%   Upwr     - {n_sensors ^ 2 x sl_rnk} SL subspace basis matrix.
%              Columns of Upwr span the VC subspace
%              If empty, compute from G2.
%   seed_ind - int scalar; index of seed source
%   cp_part  - string; can be either 'real', 'imag' or 'full'
%              determines which part of CT should be used for
%              signal subspace extraction
%   is_fast  - boolean; if true, use fast but inaccurate PSIICOS_ScanFast
%              procedure for measuring subspace correlations;
%              Inaccuracy appears because we use unprojected from SL
%              forward model in the fast procedure;
%              default = true;
% OUTPUTS:
%   corr.data  - {1 x n_sources * (n_sources - 1) / 2} vector of
%                correlations between topographies and
%                signal subspace
%   corr.IND   - {n_sources * (n_sources - 1) / 2 x 2} matrix; index mapping
%                between upwer triangle and linear indexing;
%                Usage:
%                  [~,p] = max(corr.data); [i,j] = IND(p) % find maxcorr pair
%   Cpvec      - projected away from the VC subspace sensor space cross-spectral
%                matrix
%   Upwr       - SL subspace projector.
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    %% Preparatory steps
    import ps.ProjectFromSlComplete

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
    end

    if nargin < 4
        sig_rnk = 0;
    end

    if(isempty(sl_rnk))
        if(isempty(Upwr))
            sl_rnk = 350;
        else
            sl_rnk = size(Upwr, 2);
        end
    end
    % ------------------- %


    n_sources = size(G2, 2) / 2; % two topography columns per each source
    n_sensors_C = sqrt(size(CT,1));
    n_sensors_G = size(G2,1);

    % --------------------- normalize forward matrix --------------------- %
    for i = 1:n_sources
        range_i = i * 2 - 1 : i *  2;
        G2(:, range_i(1)) = G2(:, range_i(1)) / norm(G2(:, range_i(1)));
        G2(:, range_i(2)) = G2(:, range_i(2)) / norm(G2(:, range_i(2)));
    end
    % -------------------------------------------------------------------- %

    % ---------------- check if matrix sizes are ok. --------------- %
    assert(  n_sensors_C == fix(n_sensors_C),...
            ['NONSQUARE NUMBER OF ROWS IN CT: size(CT) = ', num2str(size(CT))] );

    assert(n_sensors_C == n_sensors_G,...
           ['INCONSISTENT NUMBER OF SENSORS IN CT AND G2: ',...
             num2str(n_sensors_C), ' vs ', num2str(n_sensors_G)]);
    % -------------------------------------------------------------- %

    % -- project cross-spectrum timeseries matrix from signal-leakage subspace -- %
    if isempty(Upwr) && sl_rnk
        [Cpvec, Upwr] = ProjectFromSlComplete(CT, G2, sl_rnk);
    elseif sl_rnk % use the existing matrix if profided
        assert(n_sensors_C ^ 2 == size(Upwr,1),...
               ['INCONSISTENT SIZES: size(CT,1) = ',...
                num2str(size(CT,1)), ', size(Upwr,1) = ', num2str(size(Upwr,1))]);
        Cpvec = CT - Upwr * (Upwr' * CT);
    elseif ~sl_rnk
        Cpvec = CT;
    end
    % --------------------------------------------------------------------------- %

    % --------------- extract signal subspace --------------- %
    if sig_rnk

        if strcmp(cp_part, 'real')
            [uc,~,~] = svd(real(Cpvec), 'econ');
        elseif strcmp(cp_part, 'imag')
            [uc,~,~] = svd(imag(Cpvec), 'econ');
        elseif strcmp(cp_part, 'full')
            [uc,~,~] = svd(Cpvec, 'econ');
        end

        Cp = uc(:, 1:sig_rnk);
    elseif sig_rnk == 0
        Cp = sum(Cpvec, 2);

        if strcmp(cp_part, 'real')
            Cp = real(Cp);
        elseif strcmp(cp_part, 'imag')
            Cp = imag(Cp);
        elseif strcmp(cp_part, 'full')
            Cp = Cp;
        end

        Cp = Cp / norm(Cp, 'fro');
    else
        fprintf('ERROR: PSIICOS: Signal space rank %f\n is not valid', sig_rnk);
        return;
    end
    % ------------------------------------------------------- %


    if is_fast
        import ps.PSIICOS_ScanFast
        [corr.data, corr.IND] = PSIICOS_ScanFast(G2, Cp);

        if ~isempty(seed_ind)
            % Add 1 as a placeholder for node coherence with itself
            seed_indices = corr.IND(:,1) == seed_ind | corr.IND(:,2) == seed_ind;
            corr.data = corr.data(seed_indices);
            corr.data = [corr.data(1:seed_ind); 0; corr.data(seed_ind + 1 : end)];
            corr.IND = corr.IND(seed_indices,:);
            corr.IND = [corr.IND(1:seed_ind, :); [seed_ind, seed_ind]; corr.IND(seed_ind + 1 : end, :)];
        end
    else
        [corr.data, corr.IND] = honest_corrs(G2, Upwr, Cp, seed_ind, cp_part);
    end

end


function [Cs, IND] = honest_corrs(G2, Upwr, Cp, seed_ind, cp_part)
% -------------------------------------- %
% Compute subspace correlations
% -------------------------------------- %

    import ups.indUpperDiag2mat
    n_src = size(G2, 2) / 2; % two topography columns per each source

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
                Cs(p) = get_ij_subcorr(G2, Upwr, Cp, i_src, j_src, cp_part);

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
                Cs(i_src) = get_ij_subcorr(G2, Upwr, Cp, i_src, seed_ind, cp_part);
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
        seed_indices = IND(:,1) == seed_ind | IND(:,2) == seed_ind;
        IND = IND(seed_indices, :);
        IND = [IND(1:seed_ind-1, :); [seed_ind, seed_ind]; IND(seed_ind:end, :)];
    end
end


function cs = get_ij_subcorr(G2, Upwr, Cp, i_src, j_src, cp_part)
% --------------------------------------------------------------- %
% Get subspace correlation for [i_src, j_src] network
% --------------------------------------------------------------- %
    import ps.subcorr

    assert(i_src ~= j_src,...
          ['ERROR: get_ij_subcorr: i_src = j_src for i_src = ', num2str(i_src),...
          'j_src = ', num2str(j_src)])

    ai = G2(:, i_src * 2 - 1 : i_src * 2);
    aj = G2(:, j_src * 2 - 1 : j_src * 2);
    Gij = kron(ai, aj);
    Gji = kron(aj, ai);

    if strcmp(cp_part,'real')
        G_re = Gij + Gji;
        G_re = G_re - Upwr * (Upwr' * G_re);
        subc = subcorr(G_re, Cp);
    elseif strcmp(cp_part, 'imag')
        G_im = Gij - Gji;
        subc = subcorr(G_im, Cp);
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
