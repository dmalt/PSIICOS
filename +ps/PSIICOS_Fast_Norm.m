function [Cs, IND, Cp, Upwr] = PSIICOS_Fast_Norm(C, G2dU, Rnk, Upwr)
% ----------------------------------------------------------------------
% PSIICOS with normalization by power
% ----------------------------------------------------------------------
% FORMAT:
%   function [Cs, IND, Cp, Upwr] = PSIICOS_Fast_Norm(C, G2dU, Rnk, Upwr)
% INPUT:
%  C          - sensor space cross-spectral matrix
%  G2dU       - forward model matrix such that each source is served
%               by two columns of this matrix corresponding to the
%               topographies of dipoles in the tangential plane
%  Rnk        - scalar; limits the dimension of Volume Conduction subspace
% OUTPUT:
%  Cs         - the upper diagonal of source space cross-spectrum
%  IND        - vector of pair of dipole pairs corresponding to the
%               elements of Cs
%  Cp         - projected away from the VC subspace sensor space
%               cross-soectral matrix
%  Upwr       - VC subspace basis matrix. Columns of Upwr span the
%               VC subspace
% _______________________________________________________________________

    import ps.ProjectorOnlyAwayFromPowerComplete
    import ps.PSIICOS_ScanFast

    if(nargin < 3)
        Rnk  = 350;
    end;

    if(nargin < 4)
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
    %    [Cpvec, Upwr] = ProjectAwayFromVC(C(:), G2dU,Rnk);
    %    [Cpvec, Upwr] = ProjectAwayFromPowerFixedOr(C(:), G2dU,Rnk);
        [Upwr] = ProjectorOnlyAwayFromPowerComplete(G2dU, Rnk);
    end;

    c = Upwr' * C(:);
    Cpvec  = C(:) - Upwr * c;

    Cvcvec  = Upwr * c;
    % Cvcvec  = C(:);
    sen_pwr = diag(C);
    % PP = sen_pwr * sen_pwr';
    % size(PP)
    

    Cp = reshape(Cpvec, size(C,1), size(C,2));
    Cvc = reshape(Cvcvec, size(C,1), size(C,2));

    %% normalize forward matrix
    for i = 1 : Nsrc
        range_i = i * 2 - 1 : i * 2;
        G2dU(:,range_i(1)) = G2dU(:, range_i(1)) / norm(G2dU(:, range_i(1)));
        G2dU(:,range_i(2)) = G2dU(:, range_i(2)) / norm(G2dU(:, range_i(2)));
    end;

    Ps = zeros(Nsrc, 1);
    for i = 1 : Nsrc
        range_i = i * 2 - 1 : i * 2;
        ai = G2dU(:, range_i)';
        cs = ai * Cvc * ai';
        [~, s, ~] = svd(cs);
        Ps(i) = (s(1,1));
    end;

    %% scan all pairs with efficient vectorised implementation
    [Cs, IND] = PSIICOS_ScanFast(G2dU, Cp(:));
    Cs = sqrt(Cs);

    PPS = Ps * Ps';
    PPSUpper = zeros(1, Nsrc * (Nsrc - 1) / 2);
    k = 1;
    for i = 1 : Nsrc
        for j = i + 1 : Nsrc
            PPSUpper(k) = sqrt(PPS(i, j));
            k = k + 1;
        end;
    end;

    norm(Cs)
    disp('PPSUpper norm:');
    norm(PPSUpper)
    Cs = Cs ./ PPSUpper;
end
