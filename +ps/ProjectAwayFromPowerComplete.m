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
    end;

    Ns = size(G2dU, 2) / 2; % two topography columns per each source of the grid
    % perform projection of the coherence matrix away from the power only
    Nch = size(G2dU, 1);
    % component

    % project for each potential source
    fprintf('Collecting power subspace of coherence span...\n');
    A = zeros(Nch ^ 2, Ns * 3);
    for i=1:Ns
         gi = G2dU(:,2 * i - 1);
         v = gi * gi';
         A(:, 3 * i - 2) = v(:) / norm(v(:));
         gj = G2dU(:, 2 * i);
         v = gj * gj';
         A(:, 3 * i - 1) = v(:) / norm(v(:));
         v = gi * gj' + gj * gi';
         A(:, 3 * i) = v(:) / norm(v(:));
    end;

    fprintf('Finding eigen space...\n');
    % AA = A * A';
    % [u s] = eigs(AA, PwrRnk);
    [u s] = svd(A, 'econ');
    ds = diag(s);
    Upwr = u(:, 1:PwrRnk);

    if(size(CT,1) ~= size(Upwr,1))
        CTpvec  = CT(:) - Upwr * (Upwr' * CT(:));
        CTp = reshape(CTpvec, size(CT));
    else
        CTp  = CT - Upwr * (Upwr' * CT);
    end;
end
