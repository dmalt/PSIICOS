% function [indep_topo, c_ss_hat, PVU, SubC, INDrap, Cp, Upwr] = T_PSIICOS(C, G2dU, Rnk, Upwr)
function [INDrap, Cp, Upwr, Cs] = T_PSIICOS(C, G2dU, Rnk, SigRnk, Upwr)
% --------------------------------------------------------------------------------------------------
% Project from VC and do thresholding on correlations of sources with the cross-spectrum
% --------------------------------------------------------------------------------------------------
% FORMAT:
%   [indep_topo, c_ss_hat, PVU, SubC, INDrap, Cp, Upwr] = T_PSIICOS(C, G2dU, RAPIts, Rnk, Upwr) 
% INPUTS:
%   C        - {N_sensors_reduced x N_sensors_reduced} sensor-space cross-spectral matrix
%   G2dU     - {N_sensors_reduced x N_sources} forward model matrix such that each source 
%              is served by two columns of this matrix corresponding to the \
%              topographies of dipoles in the tangential plane
%   RAPIts   - scalar; number of iterations for the algorithm to perform
%   Rnk      - scalar; rank of Volume Conduction subspace. The bigger this value
%              is the more data will be removed by the projection from VC. On the
%              contrary, the smaller it is the more VC-related activity will remain in the data.
%   Upwr     - {N_sensors_reduced ^ 2 x Rnk} VC subspace basis matrix. 
%              Columns of Upwr span the VC subspace
% OUTPUTS:
%   indep_topo - {N_sources ^ 2 x 2 * RAPIts} matrix of oriented topographies of connected
%                pairs recovered by the algorithm
%   c_ss_hat   - {rap x 2} 
%   PVU        - {1 x RAPIts}; array of percentages of variance unexplained for each
%                algorithm iteration.
%   SubC       - 
%   Cp         - projected away from the VC subspace sensor space cross-spectral
%                matrix
%   Upwr       - VC subspace basis matrix. Columns of Upwr span the VC subspace
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    %% Preparatory steps
    if(nargin < 4)
        Rnk  = 350;
    end;

    if(nargin < 5)
        Upwr  = [];
    end;

    if(isempty(Rnk))
        if(isempty(Upwr))
            Rnk = 350;
        else
            Rnk = size(Upwr, 2);
        end;
    end;

    Nsrc = size(G2dU, 2) / 2; % two topography columns per each source of the grid
    Nch = size(G2dU, 1);

    %% perform projection of the coherence matrix away from the power only
    if(isempty(Upwr))
        % [Cpvec, Upwr] = ProjectAwayFromPowerFixedOr(C(:), G2dU, Rnk);
        [Cpvec, Upwr] = ProjectAwayFromPowerComplete(C, G2dU, Rnk);
    else % use the existing matrix if profided
        c = Upwr' * C;
        Cpvec  = C - Upwr * c;
    end;

    %k =1;
    % for rnk = 1:10:size(Upwr,2)
    %     Cimp  = imag(C(:))-Upwr(:,1:rnk)*(Upwr(:,1:rnk)'*imag(C(:)));
    %     nrmim(k) = norm(Cimp(:))/norm(imag(C(:)));
    %     k = k+1;
    %end;

    % Cp = reshape(Cpvec, size(C, 1), size(C, 2));
    [uc,sc,vc] = svd(Cpvec, 'econ');
    % SigRnk = 2;
    Cp = uc(:,1:SigRnk);
    % Cp = sum(Cpvec,2);
    % Cp = Cp / norm(Cp, 'fro');
    % Cp = Cpvec;
    %% normalize forward matrix
     
     for i = 1:Nsrc
         range_i = i * 2 - 1 : i *  2;
         G2dU(:, range_i(1)) = G2dU(:, range_i(1)) / norm(G2dU(:, range_i(1)));
         G2dU(:, range_i(2)) = G2dU(:, range_i(2)) / norm(G2dU(:, range_i(2)));
     end;

    %% scan all pairs with efficient vectorised implementation
    Urap = [];
    range2 = 1:2;
    indep_topo = zeros(prod(size(Cp)), 2);
    c_ss_hat = zeros(1, 2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % for rap = 1:RAPIts
        % Look at the topography of a pair that is
        % most correlated with the cross-spectrum
    % Cp = Cp / norm(Cp, 'fro');
    corr_threshold = 0.44;
    % [Cs(1,:), IND, Cs0] = PSIICOS_ScanFast(G2dU, Cp);
    Cs = zeros(1, Nsrc * (Nsrc - 1) / 2);
    p = 1;
    for iSrc = 1:Nsrc
        for jSrc = iSrc + 1 : Nsrc
            ai = G2dU(:,iSrc * 2 - 1 : iSrc * 2);
            aj = G2dU(:,jSrc * 2 - 1 : jSrc * 2);
            Gij = kron(ai,aj);
            Gji = kron(aj,ai);
            G = [Gij, Gji];
            subc = subcorr(G, Cp);
            Cs(1, p) = subc(1);
            p = p + 1;
        end
    end
    IND = UpperDiagToPairs(Nsrc);
    % [val_max ind_max] = max(Cs(rap,:));
    ind_threshold = find(sqrt(Cs * 2) > corr_threshold);
    pair_max = IND(ind_threshold,:);
    i = IND(ind_threshold, 1); 
    j = IND(ind_threshold, 2);
    INDrap(:,1) = i;
    INDrap(:,2) = j;
    % range_i = i * 2 - 1 : i * 2;
    % range_j = j * 2 - 1 : j * 2;
    % ai = G2dU(:, i * 2 - 1 : i * 2);
    % aj = G2dU(:, j * 2 - 1 : j * 2);
    % cs = ai' * Cprap * aj;
    % % [u s v] = svd(cs);
    % % % u and v are complex but the orientation vectors of dipoles are physical and
    % % % therefore can not be defined over the field of complex numbers. 
    % % % Do this trick to force SVD to real 2x1 vectors(first left, ten right) to find 
    % % % the real orientations

    % % --- Do we realy need this???? --- %
    % csr = real(cs);
    % csi = imag(cs);
    % [uL sL vL] = svd([csr csi]);
    % [uR sR vR] = svd([csr; csi]);
    % ai_or = ai * uL(:,1);
    % aj_or = aj * vR(:,1);
    % % --- This is correct but very slow ------ %
    % % [uL, vR, f] = FindOr(Cprap, ai, aj);
    % % ai_or = ai * uL(:,1);
    % % aj_or = aj * vR(:,1);
    % % ---------------------------------------- %
    
    % % -------------------------------- %
    % % note that norm(ai_or '* Cprap * aj_or) = s(1,1)= s_max and therefore the
    % % fast scan implemented is valid
    % qij = ai_or * aj_or';
    % qji = aj_or * ai_or';
    % qijp = qij(:) - Upwr * (Upwr' * qij(:));
    % qjip = qji(:) - Upwr * (Upwr' * qji(:));
    % qp = [qijp, qjip];
    % c_ss_hat(1, range2) = (pinv(qp) * Cprap(:))';
    % SubC(rap) = subcorr(Cprap(:), qp);
    % % Project the cross-spectrum away from the most 
    % % correlated source
    % Cprap_vec = Cprap(:) - qp * (pinv(qp) * Cprap(:));
    % Cprap = reshape(Cprap_vec, size(Cprap));
    % % Recovered independent topographies:
    % indep_topo(:,range2) = qp;
    % %%%%%%%%%%%%%%%%%%%%%%%%%
    % PVU(rap) = norm(Cp(:) - indep_topo * c_ss_hat') / norm(Cp(:));
    % % Indices of connected sited on a source level:
    % INDrap(rap, 1) = i; 
    % INDrap(rap, 2) = j; 
    % %%%%%%%%%%%%%%%%%%%
    % range2 = range2 + 2;
  % end

    

