function [U, S] = ComputeSlSvd(G, normalize, loose)
% --------------------------------------------------------------------------- %
% Compute SVD of the signal leakage matrix
% --------------------------------------------------------------------------- %
% FORMAT:
%   [U, S] = ComputeSlSvd(G, normalize, loose)
% INPUT:
%   G      - {n_sensors x n_sources * 2} matrix;
%             MEG loose orientation forward operator
%             projected to tangential plane
%   normalize - boolean;
%               if true, normalize columns before computing SVD
%   loose     - boolean;
%               if false, treat each of two dipoles per source location
%               as a separate source
% OUTPUT:
%   U      - {n_sensors ^ 2 x n_sensors ^ 2}
%            matrix of left singular vectors of signal leakage matrix
%   S      - {n_sensors ^ 2 x n_sensors ^ 2}
%            diagonal matrix of singular values of the signal leakage matrix
% ___________________________________________________________________________ %


    A = +ps.assemble_sl_matrix(G, normalize, loose);
    fprintf('Finding eigen space...\n');
    % A = [A, ones(size(A, 1), 1) * norm(A)];

    n_sen = size(G, 1);
    n_symm = n_sen * (n_sen + 1) / 2;
    n_antisymm = n_sen * (n_sen - 1) / 2;
    % set indices for lower triang part in vectorized matrix
    A_tril = convert2tril(A);
    [U_tril, S_tril] = svd(A_tril);

    U_symm = convert2mat(U_tril, 'symm');
    U_antisymm = convert2mat(U_tril(:, 1:n_antisymm), 'antisymm');
    U = [U_symm, U_antisymm];

    S = zeros(size(A));
    S(1:n_symm, :) = S_tril;

    fprintf('Done\n');
end
