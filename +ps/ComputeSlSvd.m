function [U, S] = ComputeSlSvd(G2, normalize, loose)
% --------------------------------------------------------------------------- %
% Compute SVD of the signal leakage matrix
% --------------------------------------------------------------------------- %
% FORMAT:
%   [U, S] = ComputeSlSvd(G2, normalize, loose)
% INPUT:
%   G2      - {n_sensors x n_sources * 2} matrix;
%             MEG loose orientation forward operator
%             projected to tangential plane
%   normalize - boolean;
%               if true, normalize columns before computing SVD
%   loose     - boolean;
%               if false, treat each of two dipoles per source location
%               as a separate source
% OUTPUT:
%   U      - {n_sensors ^ 2 x rank(SL matrix)}
%            matrix of left singular vectors of signal leakage matrix
%   S      - {rank{SL matrix} x rank(SL matrix)}
%            diagonal matrix of singular values of the signal leakage matrix
% ___________________________________________________________________________ %

    n_src = size(G2, 2) / 2; % in G2 we have two topography columns per source
    n_sen = size(G2, 1);

    % project for each potential source
    fprintf('Collecting signal leakage subspace...\n');
    if loose
        A = zeros(n_sen ^ 2, n_src * 3);
        for i = 1:n_src
             gi = G2(:, 2 * i - 1);
             gj = G2(:, 2 * i);

             A(:, 3 * i - 2) = kron(gi, gi);
             A(:, 3 * i - 1) = kron(gj, gj);
             A(:, 3 * i) = kron(gi, gj) + kron(gj, gi);
         end
    else
        A = zeros(n_sen ^ 2, n_src * 2);
        for i = 1 : n_src * 2
             gi = G2(:, i);
             A(:, i) = kron(gi, gi);
         end
    end

     if normalize
         A = A ./ sum(A .* 2, 1); % normalize each column
     end

    fprintf('Finding eigen space...\n');
    [U, S] = svd(A, 'econ');
end
