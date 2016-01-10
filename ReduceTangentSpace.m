function Gain_reduced = ReduceTangentSpace(Nsites, Gain, ChUsed)
% ----------------------------------------------------------------
% ReduceTangentSpace: takes takes as input Gain matrix with three 
% dipoles for each source location and reduces them to two
% laying in a tangent plane
% ----------------------------------------------------------------
% FORMAT:
%   Gain_reduced = ReduceTangentSpace(Nsites, Gain, ChUsed) 
% INPUTS:
%   Nsites        - number of source locations on cortex
%   Gain          - {Nchannels x 3 * Nsources} forward operator 
%                   gain matrix
%   ChUsed        - {1 x NumOfChannelsUsed} array of indeces of
%                     channels that are left for analysis; 
%                     to use gradiometers only go for
%                     ChUsed = 1:306; ChUsed(3:3:end) = []; 
% OUTPUTS:
%   Gain_reduced  - {Nchannels x 2 * Nsources} reduced forward 
%                     operator matrix with only two dipols per
%                     source location confined to a tangent plane
% ________________________________________________________________
% Dmitrii Altukhov, dm.altukhov@ya.ru

    range = 1:2;
    for i=1:Nsites
        g = [Gain(ChUsed, 1 + 3 * (i - 1)) ...
             Gain(ChUsed, 2 + 3 * (i - 1)) ...
             Gain(ChUsed, 3 + 3 * (i - 1))];
        [u sv v] = svd(g);
        gt = g * v(:, 1:2);
        Gain_reduced(:,range) = gt * diag(1 ./ sqrt(sum(gt .^ 2, 1)));                                                                                                                                                                                                                              
        range = range + 2;
    end;
