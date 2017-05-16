function [INDrap, Cpvec, Upwr, corr, IND] = E_PSIICOS(trials, G2dU, rel_threshold, Rnk, SigRnk, Upwr)
% --------------------------------------------------------------------------------------------------
% PSIICOS for envelopes correlation
% --------------------------------------------------------------------------------------------------
% FORMAT:
%   [INDrap, Cp, Upwr, corr] = E_PSIICOS(C, G2dU, corr_threshold, Rnk, SigRnk, Upwr) 
% INPUTS:
%   C        - {N_sensors_reduced * N_sensors_reduced x n_Times} 
%              sensor-space cross-spectral matrix
%   G2dU     - {N_sensors_reduced x N_sources} forward model matrix 
%              such that each source is served by two columns
%              of this matrix corresponding to the topographies of dipoles
%              in the tangential plane
%   Rnk      - scalar; rank of signal leakage subspace. The bigger this value
%              the more data will be removed by the projection from SL. On the
%              contrary, the smaller it is the more SL-related activity will
%              remain in the data. Recommended values are between 350 and 500
%   SigRnk   - scalar; number of components left after dimensionality reduction
%              of signal space
%              if SigRnk = 0, use mean of cross-spectrum across the time domain 
%   Upwr     - {N_sensors_reduced ^ 2 x Rnk} VC subspace basis matrix. 
%              Columns of Upwr span the VC subspace
% OUTPUTS:
%   INDrap     - {nActivePairs x 2} matrix of upper-diagonal indices of pairs 
%                with coherence levels above corr_threshold
%   Cp         - projected away from the VC subspace sensor space cross-spectral
%                matrix
%   Upwr       - VC subspace basis matrix. Columns of Upwr span the VC subspace
%   corr       - {1 x Nsrc * (Nsrc - 1) / 2} vector of correlation between
%                topographies and signal subspace
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

% Outline of the algorithm:
% 1) for each trial compute C
% 2) project C from leakage
% 3) for each C estimate lambda = Ai*Aj. Thereby we get products for each trial and each time slice
% 4) estimate powers using DICS (We get Ai for each source)
% 5) compute envelopes correlation averaging across time and maybe trials
    import ps.ProjectAwayFromPowerComplete
    import ps.PSIICOS_ScanFast
    import ups.threshold_connections
    import ups.CrossSpectralTimeseries

    %% Preparatory steps
    % if(nargin < 5)
    %     Rnk  = 350;
    % end;

    if(nargin < 6)
        Upwr  = [];
    end;

    if nargin < 5
        SigRnk = 0;
    end

    if nargin < 4
        if(isempty(Upwr))
            Rnk = 350;
        else
            Rnk = size(Upwr, 2);
        end;
    end;

    Nsrc = size(G2dU, 2) / 2; % two topography columns per each source of the grid
    % Nch = size(G2dU, 1);
    % perform projection of the coherence matrix away from the power only
    % n_sensors_C = sqrt(size(C,1));
    % n_sensors_G = size(G2dU,1);
    % assert(  n_sensors_C == fix(n_sensors_C),...
    %         ['NONSQUARE NUMBER OF ROWS IN C: size(C) = ', num2str(size(C))] ); 

    % assert(n_sensors_C == n_sensors_G,...
    %        ['INCONSISTENT NUMBER OF SENSORS IN C AND G2dU: ',...
    %          num2str(n_sensors_C), ' vs ', num2str(n_sensors_G)]);

    [n_sensors, n_times, n_trials] = size(trials);

    % ------------------- hilbert-transform trials ----------------- %
    Xfft = fft(trials, [], 2);
    h  = zeros(1, n_times); % nx1 for nonempty. 0x0 for empty.
    if n_times > 0 && 2 * fix(n_times / 2) == n_times
    % even and nonempty
    h([1 n_times / 2 + 1]) = 1;
    h(2:n_times / 2) = 2;
    elseif n_times > 0
    % odd and nonempty
    h(1) = 1;
    h(2:(n_times + 1) / 2) = 2;
    end
    HF = repmat(h, [n_sensors, 1, n_trials]);
    XH = ifft(Xfft .* HF, [], 2);
    Xph = XH; %./(abs(XH)+0.0001*mean(abs(XH(:))));
    XphConj = conj(Xph);
    % -------------------------------------------------------------- %
    clear XH;


    k = 1;
    C_tr = zeros(n_sensors ^ 2, n_times, n_trials);
    for i_trial = 1:n_trials
        for i_ch = 1:n_sensors
            C_tr(k:k + n_sensors - 1,:,:) = bsxfun(@times,...
                                               Xph(1:n_sensors,:,:),...
                                               XphConj(i_ch,:,:));
            k = k + n_sensors;
        end
    end
    clear XphConj;
    
    for i = 1:n_times
        for j = 1:n_trials
            [C_tr(:,i,j), Upwr, ~] = ProjectAwayFromPowerComplete(C_tr(:,i,j),...
                                                                  G2dU, Rnk, Upwr);
        end
    end

    % Cp = Cp / norm(Cp, 'fro');
    C_tr = C_tr * 1e23;

    %% normalize forward matrix
     
     for i = 1:Nsrc
         range_i = i * 2 - 1 : i *  2;
         G2dU(:, range_i(1)) = G2dU(:, range_i(1)) / norm(G2dU(:, range_i(1)));
         G2dU(:, range_i(2)) = G2dU(:, range_i(2)) / norm(G2dU(:, range_i(2)));
     end;

    AiAj = zeros(1, (n_sensors ^ 2 - n_sensors) / 2);
    for k = 1:n_times
        for l = 1:n_trials
            AiAj(:) = AiAj + PSIICOS_ScanFast(G2dU, C_tr(:,k,l));
        end
    end

    AiAj = AiAj / (n_times * n_trials);
    % Etimate amplutide means Ai and Aj here
     
    CT = ups.CrossSpectralTimeseries(trials, false);
    for i = 1:n_times
        CT(:,i) =  Upwr * (Upwr' * C(:,i));
    end


    Ps = zeros(Nsrc, n_times);

    A = G2dU';

    for t = 1:n_times
        C = reshape(CT(:,t), n_sensors, n_sensors);
        for i = 1:Nsrc
            range_i = i * 2 - 1 : i * 2;
            ai = A(range_i,:);
            cs = ai * C * ai';
            [~, s, ~] = svd(cs);
            Ps(i,t) = sqrt(s(1,1));
        end
    end;

end
