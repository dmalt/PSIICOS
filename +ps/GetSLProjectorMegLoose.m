function [Upwr, ds] = GetSLProjectorMegLoose(G2, pwr_rnk)
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

    n_src = size(G2, 2) / 2; % in G2 we have two topography columns per source
    n_sen = size(G2, 1);

    % project for each potential source
    fprintf('Collecting signal leakage subspace...\n');
    A = zeros(n_sen ^ 2, n_src * 3);
    for i = 1:n_src
         gi = G2(:,2 * i - 1);
         v = gi * gi';
         A(:, 3 * i - 2) = v(:) / norm(v(:));
         gj = G2(:, 2 * i);
         v = gj * gj';
         A(:, 3 * i - 1) = v(:) / norm(v(:));
         v = gi * gj' + gj * gi';
         A(:, 3 * i) = v(:) / norm(v(:));
     end

    fprintf('Finding eigen space...\n');
    [u, s] = svd(A, 'econ');
    ds = diag(s);
    Upwr = u(:, 1:pwr_rnk);
end
