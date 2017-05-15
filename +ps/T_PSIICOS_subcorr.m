function [INDrap, IND, Cp, Upwr, Cs] = T_PSIICOS_subcorr(C, G2dU, SL_rnk,...
                                                         sig_rnk, Upwr,...
                                                         seed_ind, cp_part)
% ---------------------------------------------------------------------------------------
% Experimential version with proper subspace correlations.
% ---------------------------------------------------------------------------------------
% FORMAT:
%   [indep_topo, c_ss_hat, PVU, SubC, INDrap, Cp, Upwr] = T_PSIICOS(C, G2dU, RAPIts, SL_rnk, Upwr) 
% INPUTS:
%   C        - {N_sensors_reduced x N_sensors_reduced} sensor-space cross-spectral matrix
%              sensor-space cross-spectral matrix
%   G2dU     - {N_sensors_reduced x N_sources} forward model matrix 
%              such that each source is served by two columns
%              of this matrix corresponding to the topographies of dipoles
%              in the tangential plane
%   SL_rnk      - scalar; rank of signal leakage subspace. The bigger this value
%              the more data will be removed by the projection from SL. On the
%   Upwr     - {N_sensors_reduced ^ 2 x SL_rnk} VC subspace basis matrix. 
%              Columns of Upwr span the VC subspace
% OUTPUTS:
%   Cp         - projected away from the VC subspace sensor space cross-spectral
%                matrix
%   Upwr       - VC subspace basis matrix. Columns of Upwr span the VC subspace
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    %% Preparatory steps
    import ps.ProjectAwayFromPowerComplete
    import ps.subcorr
    import ups.indUpperDiag2mat


    if(nargin < 5)
        Upwr  = [];
    end;

    if nargin < 4
        sig_rnk = 0;
    end

    if(isempty(SL_rnk))
        if(isempty(Upwr))
            SL_rnk = 350;
        else
            SL_rnk = size(Upwr, 2);
        end;
    end;

    Nsrc = size(G2dU, 2) / 2; % two topography columns per each source of the grid
    n_sensors_C = sqrt(size(C,1));
    n_sensors_G = size(G2dU,1);
    assert(  n_sensors_C == fix(n_sensors_C),...
            ['NONSQUARE NUMBER OF ROWS IN C: size(C) = ', num2str(size(C))] ); 

    assert(n_sensors_C == n_sensors_G,...
           ['INCONSISTENT NUMBER OF SENSORS IN C AND G2dU: ',...
             num2str(n_sensors_C), ' vs ', num2str(n_sensors_G)]);

    %% perform projection of the coherence matrix away from the power only
    if(isempty(Upwr))
        [Cpvec, Upwr] = ProjectAwayFromPowerComplete(C, G2dU, SL_rnk);
    else % use the existing matrix if profided
        assert(n_sensors_C ^ 2 == size(Upwr,1),...
               ['INCONSISTENT SIZES: size(C,1) = ',...
                num2str(size(C,1)), ', size(Upwr,1) = ', num2str(size(Upwr,1))]);
        c = Upwr' * C;
        Cpvec  = C - Upwr * c;
    end;

    % Cpvec = C;

    if sig_rnk
        [uc,~,~] = svd(Cpvec, 'econ');
        Cp = uc(:,1:sig_rnk);
    elseif sig_rnk == 0
        Cp = sum(Cpvec, 2);
        Cp = Cp / norm(Cp, 'fro');
    else
        fprintf('ERROR: PSIICOS: Signal space rank %f\n is not valid', sig_rnk);
        return;
    end
    %% normalize forward matrix
     
     for i = 1:Nsrc
         range_i = i * 2 - 1 : i *  2;
         G2dU(:, range_i(1)) = G2dU(:, range_i(1)) / norm(G2dU(:, range_i(1)));
         G2dU(:, range_i(2)) = G2dU(:, range_i(2)) / norm(G2dU(:, range_i(2)));
     end;

    if isempty(seed_ind)
        Cs = zeros(Nsrc * (Nsrc - 1) / 2,1);
    else
        Cs = zeros(Nsrc, 1);
        aj = G2dU(:, seed_ind * 2 - 1 : seed_ind * 2);
    end

    p = 1;
    for iSrc = 1:Nsrc
        if isempty(seed_ind)
            for jSrc = iSrc + 1 : Nsrc
                ai = G2dU(:,iSrc * 2 - 1 : iSrc * 2);
                aj = G2dU(:,jSrc * 2 - 1 : jSrc * 2);
                Gij = kron(ai,aj);
                Gji = kron(aj,ai);
                G = [Gij, Gji];
                subc = subcorr(G, Cp);
                Cs(p,1) = subc(1);
                p = p + 1;
            end
        else
            ai = G2dU(:,iSrc * 2 - 1 : iSrc * 2);
            Gij = kron(ai,aj);
            Gji = kron(aj,ai);
            % G = Gij + Gji;
            if strcmp(cp_part,'real')
                G_re = Gij + Gji;
                G_re = G_re - Upwr * Upwr' * G_re;
                subc = subcorr(G_re, real(Cp));
            elseif strcmp(cp_part, 'imag')
                G_im = Gij - Gji;
                G_im = G_im - Upwr * Upwr' * G_im;
                if(seed_ind == iSrc)
                    subc = rand;
                else
                    subc = subcorr(G_im, imag(Cp));
                end
            elseif strcmp(cp_part, 'full')
                G_re = Gij + Gji;
                G_re = G_re - Upwr * Upwr' * G_re;
                subc_re = subcorr(G_re, real(Cp));
                G_im = Gij - Gji;
                % G_im = G_im - Upwr * Upwr' * G_im;
                if(seed_ind == iSrc)
                    subc_im = rand;
                else
                    subc_im = subcorr(G_im, imag(Cp));
                end
                subc = sqrt(subc_re .^ 2 + subc_im .^ 2);
            else
                error(['Unknown option ', cp_part]);
            end

            % subc = subcorr(Gij, Cp);
            Cs(iSrc) = subc(1);
            disp(iSrc);
        end
    end

    IND = indUpperDiag2mat(Nsrc);
    % [val_max ind_max] = max(Cs(rap,:));
end
