function [Cs, IND Cs0] = PSIICOS_ScanFast(G2dU, Cp)
% PSIICOS_ScanFast: perform scanning algorithm to find strongest connections
%                   using projection from VC
% FORMAT:
%   [Cs, IND Cs0] = PSIICOS_ScanFast(G2dU, Cp) 
% INPUTS:
%   G2dU       - {Nsensors x Nsources} matrix of forward model
%   Cp         - {Nsensors x Time} matrix of cross-spectrum on sensors
% OUTPUTS:
%   Cs
%   IND
%   Cs0
% ___________________________________________________________________________
% Alex Ossadtchii, ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    [Nsns, Nsrc2] = size(G2dU);
    Nsrc = Nsrc2 / 2;

    if(size(Cp,1)~=Nsns)
        disp('Incomptible dimensions G2dU vs Cp');
        return;
    end

    T = zeros(1, Nsrc * (Nsrc - 1) / 2);
    D = zeros(1, Nsrc * (Nsrc - 1) / 2);
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
    %below is the optimized implementation of this:
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
    tic
    p = 1;
    ai = zeros(2, Nsns);
    for i=1:Nsrc
        % --- Take i-th location topographies ---- %
        ai = G2dU(:, i * 2 - 1 : i * 2)';       
        tmp = ai * Cp;
        cslong =  tmp * G2dU;
        cs2long = cslong .* conj(cslong);
        cs2longd = cslong(1,:) .* conj(cslong(2,:));
        cs2_11_22 = [sum(reshape(cs2long(1,:), 2, Nsrc), 1);...
                    sum(reshape(cs2long(2,:), 2, Nsrc), 1)];
        cs2_12_21 = sum(reshape(cs2longd, 2, Nsrc), 1);
        Ti = sum(cs2_11_22, 1);
        Di = prod(cs2_11_22, 1) - cs2_12_21 .* conj(cs2_12_21);
        T(p : p + Nsrc - 1 - i) = Ti(i + 1 : Nsrc);
        D(p : p + Nsrc - 1 - i) = Di(i + 1 : Nsrc);
        IND(p : p + Nsrc - i - 1, 2) = (i + 1 : Nsrc)';
        IND(p : p + Nsrc - i - 1, 1) = i;
        p = p + Nsrc - i;
    end;
    Cs = 0.5 * T + sqrt(0.25 * T .* (T) - D); 
           
    toc



