function corr = RAP_PSIICOS(CT, G2, sl_rnk, sig_rnk,...
                     Upwr, seed_ind, cp_part, is_fast, n_rap)
% --------------------------------------------------
% Recursively-applied PSIICOS
% --------------------------------------------------
% FORMAT:
%   [corr, Cpvec, Upwr] = RAP_PSIICOS(CT, G2, sl_rnk, sig_rnk, Upwr,...
    %                                 seed_ind, cp_part, is_fast, n_rap)
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
%   n_rap    - number of PSIICOS iterations.
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
%   corr       - {n_rap x 1} cell array; each cell contains PSIICOS solution for
%                the corresponding iteration (see OUTPUT section of PSIICOS docstring)
% _____________________________________________________________________________________
% Dimtrii Altukhov, dm.altukhov@ya.ru
%
    import ps.PSIICOS

for i_rap = 1:n_rap
     [corr{i_rap}, CT, Upwr] = PSIICOS(CT, G2, sl_rnk, sig_rnk, Upwr,...
                                       seed_ind, cp_part, is_fast);

     if i_rap == 1
         Upwr_ini = Upwr;
         sl_rnk = 0;
     end

     % Find strongest connection and project CT away from it
     [~, indmax] = max(corr{i_rap}.data);
     ij_max = corr{i_rap}.IND(indmax,:);

     g_i_xy = G2(:, ij_max(1) * 2 - 1 : ij_max(1) * 2);
     g_j_xy = G2(:, ij_max(2) * 2 - 1 : ij_max(2) * 2);

     G_ij = kron(g_i_xy, g_j_xy);
     G_ji = kron(g_j_xy, g_i_xy);

     GG = [G_ij, G_ji];
     GG = GG - Upwr_ini * (Upwr_ini' * GG);
     CT = CT - GG * (pinv(GG) * CT);
     % Upwr = zeros(size(Upwr)); % don't do projection on 2-nd and subseq. iters
end
