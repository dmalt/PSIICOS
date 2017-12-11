function [CTp, Upwr] = ProjectFromSlComplete(CT, G2dU, pwr_rnk)
% -----------------------------------------------------------------------
% Project vectorized cross-spectrum away from power subspace
% -----------------------------------------------------------------------
% FORMAT:
%   [CTp, Upwr] = ProjectFromSlComplete(CT, G2dU, pwr_rnk)
% INPUTS:
%   CT        - {n_sensors ^ 2 x n_times} matrix;
%               vectorized cross-spectrum timeseries
%   G2dU      - {n_sensors x n_sources} matrix;
%               forward operator
%   pwr_rnk    - scalar; rank of projector
% OUTPUTS:
%   CTp       - {n_sensors ^ 2 x n_times} matrix;
%               vectorized CP timeseries proj. from VC
%   Upwr      - {n_sensors ^ 2 x pwr_rnk} matrix;
%               projeco
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    if(nargin < 3)
        pwr_rnk = 350;
    end

    Upwr = ps.GetSLProjectorMegLoose(G2dU, pwr_rnk);

    if(size(CT,1) ~= size(Upwr,1))
        CTpvec  = CT(:) - Upwr * (Upwr' * CT(:));
        CTp = reshape(CTpvec, size(CT));
    else
        CTp  = CT - Upwr * (Upwr' * CT);
    end
end
