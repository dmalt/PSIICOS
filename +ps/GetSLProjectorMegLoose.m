function [Upwr, ds] = GetSLProjectorMegLoose(G2, pwr_rnk, normalize)
% ------------------------------------------------------- %
% Compute projector to signal leakage subspace using
% complete set of potential sources
% ------------------------------------------------------- %
% FORMAT:
%   [Upwr,ds] = GetSLProjectorMegLoose(G2, pwr_rnk)
% INPUTS:
%   G2      - {n_sensors x n_sources * 2} matrix;
%             MEG loose orientation forward operator
%             projected to tangential plane
%   pwr_rnk - scalar; rank of projector
% OUTPUTS:
%   Upwr    - {n_sensors ^ 2 x pwr_rnk} matrix
%   ds      - vector; singular values of signal leakage
%             matrix
% _______________________________________________________ %
% Dmitrii Altukhov, dm.altukhov@ya.ru
    if nargin < 3
        normalize = true;
    end
    loose = true;
    [u, s] = ps.ComputeSlSvd(G2, normalize, loose);
    ds = diag(s);
    Upwr = u(:, 1:pwr_rnk);
end
