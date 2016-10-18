function [CTp, Upwr, ds, nrmre, nrmim, nrmvc] = ProjectAwayFromPowerFixedOr(CT, G2dU, PwrRnk)
% --------------------------------------------------------------------------------------------
% ProjectAwayFromPower: 
% --------------------------------------------------------------------------------------------
% FORMAT:
%   [CTp, Upwr, ds, nrmre, nrmim, nrmvc] = ProjectAwayFromPowerFixedOr(CT, G2dU, PwrRnk) 
% INPUTS:
%   CT        - {N_sensros_reduced ^ 2 x 1} vectorized cross-spectrum
%   G2dU      - {N_sensors_reduced x N_sources} forward model matrix
%   PwrRnk    - scalar; rank of a projector null-space. default = 350
% OUTPUTS:
%   CTp       - {N_sensors_reduced ^ 2 x 1} vectorized cross-spectrum projected away 
%               from volume conduction
%   Upwr      - {N_sensors_reduced ^ 2 x PwrRnk} matrix; its columns are first PwrRnk
%               singular vectors of the matrix of power topographies in squared space
%   ds        - {PwrRnk x 1} vector of first PwrRnk singular values of power subspace
%   nrmre
%   nrmim
%   nrmvc
% _________________________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru


if(nargin < 3)
    PwrRnk = 350;
end;
if(nargin < 2)
    disp('use error\n');
    CTp = [];
    Upwr = [];
    return;
end;

Ns = size(G2dU, 2) / 2; % two topography columns per each source of the grid
% perform projection of the coherence matrix away from the power only
Nch = size(G2dU, 1);
% component

% project for each potential source
fprintf('Collecting power subspace of coherence span...\n');
A = zeros(Nch ^ 2, Ns * 2);
for i=1:Ns
     gi = G2dU(:, 2 * i - 1);
     v = gi * gi';
     A(:,2 * i - 1) = v(:) / norm(v(:));
     gj = G2dU(:, 2 * i);
     v = gj * gj';
     A(:, 2 * i) = v(:) / norm(v(:));
end;

fprintf('Finding eigen space...\n');
AA = A * A';
[u s] = eigs(AA, PwrRnk);
ds = diag(s);
Upwr = u(:, 1:PwrRnk);
if(size(CT, 1) ~= size(Upwr, 1))
    CTpvec  = CT(:) - Upwr * (Upwr' * CT(:));
    CTp = reshape(CTpvec, size(CT));
else
    CTp  = CT - Upwr * (Upwr' * CT);
end;


nrmin = [];
nrmre = [];
nrmvc = [];

if(nargout == 6)
     k = 1;
     for rnk = 1:10:size(u,2)
        CTImp  = imag(CT) - u(:,1:rnk) * (u(:,1:rnk)' * imag(CT));
        CTRep  = real(CT) - u(:,1:rnk) * (u(:,1:rnk)' * real(CT));
        CTVCp  = A - u(:,1:rnk) * (u(:,1:rnk)' * A);
        nrmim(k) = norm(CTImp(:)) / norm(imag(CT(:)));
        nrmre(k) = norm(CTRep(:)) / norm(real(CT(:)));
        nrmvc(k) = norm(CTVCp(:)) / norm(A(:));
        k = k + 1;
    end;
end;
