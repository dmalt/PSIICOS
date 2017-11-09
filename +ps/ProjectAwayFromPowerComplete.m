function [CTp, Upwr, ds] = ProjectAwayFromPowerComplete(CT, G2dU, PwrRnk)
% -----------------------------------------------------------------------
% Project vectorized cross-spectrum away from power subspace
% -----------------------------------------------------------------------
% FORMAT:
%   [CTp, Upwr, ds, nrmre, nrmim, nrmvc] = ProjectAwayFromPowerComplete(CT, G2dU, PwrRnk) 
% INPUTS:
%   CT        - {n_sensors ^ 2 x n_times} matrix;
%               vectorized cross-spectrum timeseries
%   G2dU      - {n_sensors x n_sources} matrix;
%               forward operator 
%   PwrRnk    - scalar; rank of projector
% OUTPUTS:
%   CTp       - {n_sensors ^ 2 x n_times} matrix;
%               vectorized CP timeseries proj. from VC 
%   Upwr      - {} matrix
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    if(nargin < 3)
        PwrRnk = 350;
    end
    Upwr  = ps.GetUpwrComplete(G2dU, PwrRnk);
 
    if(size(CT,1) ~= size(Upwr,1))
        CTpvec  = CT(:) - Upwr * (Upwr' * CT(:));
        CTp = reshape(CTpvec, size(CT));
    else
        CTp  = CT - Upwr * (Upwr' * CT);
    end
end
