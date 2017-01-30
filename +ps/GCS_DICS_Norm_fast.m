function [Cs, Ps, IND] = GCS_DICS_Norm_fast(C,G)
% -------------------------------------------------------
% Apply geometric correction scheme + DICS beamformer
% -------------------------------------------------------
% FORMAT:
%   [Cs, Ps, IND] = GCS_DICS(C,G) 
% INPUTS:
%   C
%   G
% OUTPUTS:
%   Cs
%   Ps
%   IND
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru
 
    Nch = size(C,1);
    iC = inv(C + 40 * trace(C) / Nch * eye(Nch));
    Ns = fix(0.5 * size(G,2)); % assume tangent space dimension of 2 

    range = 1:2;
    A = zeros(size(G'));

    for i=1:Ns
        L = G(:,range);
        A(range,:) = inv(L' * iC * L) * L' * iC;
        range = range + 2;
    end

    for i=1:Ns
        range_i = i * 2 - 1 : i * 2;
        ai = A(range_i,:);
        cs = ai * C * ai';
        [~, s, ~] = svd(cs);
        Ps(i) = (s(1,1));
    end;

    p = 1;
    Cs  = zeros(Ns * (Ns - 1) / 2, 1);
    IND = zeros(Ns * (Ns - 1) / 2, 2);
    fprintf('iDICS searching the grid... \n');
    fprintf('Reference index (Max %d) : ', Ns); 

    [Cs,IND] = GCS_ScanFast(C,G,A);
    fprintf('\n Done\n');

    Cs = sqrt(Cs);
    PPS = Ps'*Ps;
    PPSUpper = zeros(1,Ns*(Ns-1)/2);
    k=1;
    for i=1:Ns
        for j = i+1:Ns
            PPSUpper(k) = sqrt(PPS(i,j));
            k = k+1;
        end;
    end;

    Cs = Cs./PPSUpper;
end

