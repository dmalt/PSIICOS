function [Cs, IND] = PSIICOS_ScanFast(G2dU, Cp, is_imag)
% -------------------------------------------------------------------------
% Perform scanning algorithm to find strongest connections using correlation
% of cross-spectrum with the forward operator
% -------------------------------------------------------------------------
% FORMAT:
%   [Cs, IND Cs0] = PSIICOS_ScanFast(G2dU, Cp) 
% INPUTS:
%   G2dU       - {n_sensors x n_sources} matrix of forward model
%   Cp         - {n_sensors ^ 2 x n_times} or
%                {n_sensors ^ 2 x n_components}
%                matrix of timeseries or left singular vectors of
%                cross-spectrum on sensors
% OUTPUTS:
%   Cs         - {(n_sources ^ 2 - n_sources) / 2} matrix of 
%                correlations between source topographies
%                and forward operator
%   IND        - {(n_sources ^ 2 - n_sources) / 2} matrix of
%                indices to build a mapping between upper
%                triangle and (i,j) matrix indexing 
%                IND(l,:) --> [i,j]
% ___________________________________________________________________________
% Alex Ossadtchii, ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    if nargin < 3
        is_imag = false;
    end

    [Nsns, Nsrc2] = size(G2dU);
    Nsrc = Nsrc2 / 2;

    % if(size(Cp,1)~=Nsns)
    %     disp('Incomptible dimensions G2dU vs Cp');
    %     return;
    % end

    n_comp = size(Cp, 2);
    T = zeros(n_comp, Nsrc * (Nsrc - 1) / 2);
    D = zeros(n_comp, Nsrc * (Nsrc - 1) / 2);
    Ti = zeros(1, Nsrc);
    Di = zeros(1, Nsrc);
    IND = zeros(Nsrc * (Nsrc - 1) / 2, 2);
    cs2 = zeros(2, 2);
    cs = zeros(2, 2);
    tmp = zeros(2, Nsrc * 2);
    ai = zeros(2, Nsns);
    aj = zeros(Nsns, 2);
    cslong = zeros(2, Nsrc * 2);
    cs2long = zeros(2, Nsrc * 2);
    cs2longd = zeros(1, Nsrc * 2);
    cs2_11_22 = zeros(2, Nsrc);
    cs2_12_21 = zeros(1, Nsrc);
    Cs0 =  zeros(1, Nsrc * (Nsrc - 1) / 2);

    % below is the optimized implementation of this:
    % Look at each pair and estimate subspace correlation
    % between cross-spectrum and topography of this pair
    % tic
    % p = 1;
    % for i=1:Nsrc
    %     range_i = i * 2 - 1 : i * 2;
    %     ai = G2dU(:,range_i)';
    %     for j=i + 1:Nsrc
    %          range_j = j * 2 - 1 : j * 2;
    %          aj = G2dU(:, range_j)';
    %         cs = ai * Cp * aj';
    %         [u s v] = svd(cs);
    %         Cs0(p) = max(diag(s)); 
    %          p = p + 1;
    %     end;
    % end;
    % toc
    %

    p = 1;
    ai = zeros(2, Nsns);
    for iSrc = 1:Nsrc
        for iComp = 1:n_comp
            Cp_sq = reshape(Cp(:,iComp), Nsns, Nsns);
            % --- Take iSrc-th location topographies ---- %
            ai = G2dU(:, iSrc * 2 - 1 : iSrc * 2)';       
            tmp = ai * Cp_sq;
            cslong = tmp * G2dU;
            if is_imag
                cslong = imag(cslong);
            end
            cs2long = cslong .* conj(cslong);
            cs2longd = cslong(1,:) .* conj(cslong(2,:));
            cs2_11_22 = [sum(reshape(cs2long(1,:), 2, Nsrc), 1);...
                         sum(reshape(cs2long(2,:), 2, Nsrc), 1)];
            cs2_12_21 = sum(reshape(cs2longd, 2, Nsrc), 1);
            Ti = sum(cs2_11_22, 1);
            Di = prod(cs2_11_22, 1) - cs2_12_21 .* conj(cs2_12_21);
            T(iComp, p : p + Nsrc - 1 - iSrc) = Ti(iSrc + 1 : Nsrc);
            D(iComp, p : p + Nsrc - 1 - iSrc) = Di(iSrc + 1 : Nsrc);
        end
        IND(p : p + Nsrc - iSrc - 1, 2) = (iSrc + 1 : Nsrc)';
        IND(p : p + Nsrc - iSrc - 1, 1) = iSrc;
        p = p + Nsrc - iSrc;
    end;
    Cs = 0.5 * T + sqrt(0.25 * T .* T - D); 
    % Cs = sum(Cs,1);  
    Cs = max(Cs,[],1);    
end
