function [INDrap, Cpvec, Upwr, corr, IND] = T_PSIICOS(C, G2dU, rel_threshold, Rnk, SigRnk, Upwr)
% --------------------------------------------------------------------------------------------------
% Project from VC and do thresholding on correlations of sources with the cross-spectrum
% --------------------------------------------------------------------------------------------------
% FORMAT:
%   [INDrap, Cp, Upwr, corr] = T_PSIICOS(C, G2dU, corr_threshold, Rnk, SigRnk, Upwr) 
% INPUTS:
%   C        - {N_sensors_reduced * N_sensors_reduced x n_Times} 
%              sensor-space cross-spectral matrix
%   G2dU     - {N_sensors_reduced x N_sources} forward model matrix 
%              such that each source is served by two columns
%              of this matrix corresponding to the topographies of dipoles
%              in the tangential plane
%   Rnk      - scalar; rank of signal leakage subspace. The bigger this value
%              the more data will be removed by the projection from SL. On the
%              contrary, the smaller it is the more SL-related activity will
%              remain in the data. Recommended values are between 350 and 500
%   SigRnk   - scalar; number of components left after dimensionality reduction
%              of signal space
%              if SigRnk = 0, use mean of cross-spectrum across the time domain 
%   Upwr     - {N_sensors_reduced ^ 2 x Rnk} VC subspace basis matrix. 
%              Columns of Upwr span the VC subspace
% OUTPUTS:
%   INDrap     - {nActivePairs x 2} matrix of upper-diagonal indices of pairs 
%                with coherence levels above corr_threshold
%   Cp         - projected away from the VC subspace sensor space cross-spectral
%                matrix
%   Upwr       - VC subspace basis matrix. Columns of Upwr span the VC subspace
%   corr       - {1 x Nsrc * (Nsrc - 1) / 2} vector of correlation between
%                topographies and signal subspace
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    import ps.ProjectAwayFromPowerComplete
    import ps.PSIICOS_ScanFast
    import ups.threshold_connections

    %% Preparatory steps
    % if(nargin < 5)
    %     Rnk  = 350;
    % end;

    if(nargin < 6)
        Upwr  = [];
    end;

    if nargin < 5
        SigRnk = 0;
    end

    if nargin < 4
        if(isempty(Upwr))
            Rnk = 350;
        else
            Rnk = size(Upwr, 2);
        end;
    end;

    Nsrc = size(G2dU, 2) / 2; % two topography columns per each source of the grid
    % Nch = size(G2dU, 1);
    % perform projection of the coherence matrix away from the power only
    n_sensors_C = sqrt(size(C,1));
    n_sensors_G = size(G2dU,1);
    assert(  n_sensors_C == fix(n_sensors_C),...
            ['NONSQUARE NUMBER OF ROWS IN C: size(C) = ', num2str(size(C))] ); 

    assert(n_sensors_C == n_sensors_G,...
           ['INCONSISTENT NUMBER OF SENSORS IN C AND G2dU: ',...
             num2str(n_sensors_C), ' vs ', num2str(n_sensors_G)]);

    if(isempty(Upwr))
        % [Cpvec, Upwr] = ProjectAwayFromPowerFixedOr(C(:), G2dU, k);
        [Cpvec, Upwr] = ProjectAwayFromPowerComplete(C, G2dU, Rnk);
    else % use the existing matrix if profided
        assert(n_sensors_C ^ 2 == size(Upwr,1),...
               ['INCONSISTENT SIZES: size(C,1) = ',...
                num2str(size(C,1)), ', size(Upwr,1) = ', num2str(size(Upwr,1))]);
        c = Upwr' * C;
        Cpvec  = C - Upwr * c;
    end;
    % Cpvec = C;
    
    if SigRnk
        [uc,~,~] = svd(Cpvec, 'econ');
        Cp = uc(:,1:SigRnk);
    elseif SigRnk == 0
        Cp = sum(Cpvec,2);
        Cp = Cp / norm(Cp, 'fro');
    else
        fprintf('ERROR: T_PSIICOS: Signal space rank %f\n is not valid', SigRnk);
        return;
    end
    %% normalize forward matrix
     
     for i = 1:Nsrc
         range_i = i * 2 - 1 : i *  2;
         G2dU(:, range_i(1)) = G2dU(:, range_i(1)) / norm(G2dU(:, range_i(1)));
         G2dU(:, range_i(2)) = G2dU(:, range_i(2)) / norm(G2dU(:, range_i(2)));
     end;

    %% scan all pairs with efficient vectorised implementation
    % Urap = [];
    % range2 = 1:2;
    % indep_topo = zeros(prod(size(Cp)), 2);
    % c_ss_hat = zeros(1, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Look at the topography of a pair that is
    % most correlated with the cross-spectrum
    [corr(1,:), IND] = PSIICOS_ScanFast(G2dU, Cp);

    % corr_min  = min(corr);
    % corr_max = max(corr);
    % corr_delta = corr_max - corr_min;
    % corr_threshold = corr_min + corr_delta * rel_threshold;
    % % ind_threshold = find(sqrt(corr * 2) > corr_threshold);
    % ind_threshold = find(corr > corr_threshold);
    % pair_max = IND(ind_threshold,:);
    % i = IND(ind_threshold, 1); 
    % j = IND(ind_threshold, 2);
    INDrap = threshold_connections(corr, rel_threshold, IND);
    % INDrap(:,1) = i;
    % INDrap(:,2) = j;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
