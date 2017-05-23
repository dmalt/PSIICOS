function [Cs, IND, Upwr] = PSIICOS_DICS(C, G, lambda, rnk, Upwr)
% -------------------------------------------------------
% PSIICOS scan with DICS inverse operator
% -------------------------------------------------------
% FORMAT:
%   [Cs, IND] = PSIICOS_DICS(C, G, lambda) 
% INPUTS:
%   C        -
%   G
%   lambda 
% OUTPUTS:
%   Cs
%   IND
% ________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru
    import ps.PSIICOS_ScanFast
    import ups.DICS
    import ps.ProjectAwayFromPowerComplete

    if nargin < 5
        Upwr = [];
    end
    if isempty(Upwr)
        [C_proj, Upwr] = ProjectAwayFromPowerComplete(C, G, rnk);
    else
        temp = Upwr' * C;
        C_proj = C - Upwr * temp;
    end

    Nch = sqrt(size(C,1));
    C_sq = reshape(C, Nch, Nch);

    A = DICS(C_sq, G, lambda);

    [Cs,IND] = PSIICOS_ScanFast(A',C_proj, false);
end
