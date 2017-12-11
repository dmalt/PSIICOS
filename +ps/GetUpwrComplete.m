function Upwr = GetUpwrComplete(G2dU, pwr_rnk)
% --------------------------------------------- %
% Compute Upwr matrix from power subspace using a
% complete set of potential sources 
% --------------------------------------------- %
% FORMAT:
% Upwr = ProjectorComplete(G2dU, PwrRnk)
% INPUTS:
%   G2dU      - {n_sensors x n_sources} matrix;
%               forward operator 
%   PwrRnk    - scalar; rank of projector
% OUTPUTS:
%   Upwr      - {n_sensors ^ 2 x pwr_rnk} matrix
% _____________________________________________ %
% Dmitrii Altukhov, dm.altukhov@ya.ru
    Ns = size(G2dU, 2) / 2; % two topography columns per each source of the grid
    % perform projection of the coherence matrix away from the power only
    Nch = size(G2dU, 1);

    % project for each potential source
    % component
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
    [u, s] = svd(A, 'econ');
    ds = diag(s);
    Upwr = u(:, 1:pwr_rnk);
end
