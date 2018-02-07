function Upwr = proj_tril(G2dU, PwrRnk, W_inv)

    r_g = rank(G2dU);
    rr = r_g * (r_g + 1) / 2;
    tril_mask = true(r_g);
    tril_mask = tril(tril_mask);
    tril_mask = tril_mask(:);

    if(nargin < 2)
        PwrRnk = 350;
    end
    if(nargin < 2)
        disp('use error\n');
        % CTp = [];
        Upwr = [];
        return;
    end

    Ns = size(G2dU, 2) / 2; % two topography columns per each source of the grid
    % perform projection of the coherence matrix away from the power only
    Nch = size(G2dU, 1);
    % component

    % project for each potential source
    fprintf('Collecting power subspace of coherence span...\n');
    A = zeros(rr, Ns * 3);
    % W_inv = inv(W);
    for i = 1 : Ns
         gi = G2dU(:, 2 * i - 1);
         v = kron(gi, gi);
         vv = v(tril_mask,:);
         % vv = W_inv * vv;
         A(:,3 * i - 2) = vv / norm(vv);
         gj = G2dU(:,2 * i);
         v = kron(gj, gj);
         vv = v(tril_mask,:);
         % vv = W_inv * vv;
         A(:,3 * i - 1) = vv / norm(vv);
         v = kron(gi, gj) + kron(gj, gi);
         vv = v(tril_mask,:);
         % vv = W_inv * vv;
         A(:, 3 * i) = vv / norm(vv);
    end
    A = W_inv * A;

    fprintf('Finding eigen space...\n');
    AA = A * A';
    [u, s] = eigs(AA, PwrRnk);
    ds = diag(s);
    Upwr = u(:, 1:PwrRnk);
end
