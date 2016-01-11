function ConData = ReduceDimensions(ConData, ChUsed, bUseHR)
% -------------------------------------------------------------------------
% ReduceDimensions: for each condition recording in ConData performs PCA
% on a gain matrix to confine three dipoles to two in a tangential plane 
% and then reduces the dimesionality of sensor space via PCA (using spm_svd).
% Effectively this means that we are left with somewhat about 40 artificial 
% sensors instead of 204 if we were using only gradiometers
% --------------------------------------------------------------------------
% FORMAT:
%   ReduceDimensions(ConData, ChUsed, bUseHR) 
% INPUTS:
%   ConData        - {1 x N_conditions_total} cell array with each element 
%                    containing a structure for a single condition 
%                    recording from protocol. Gain matrices are stored in
%                    ConData{i}.HM_LR and ConData{i}.HM_HR 
%   ChUsed         - {1 x NumOfChannelsUsed} array of indeces of
%                    channels that are left for analysis; 
%                    to use gradiometers only go for
%                    ChUsed = 1:306; ChUsed(3:3:end) = [];
%   bUseHR         - boolean flag; if True high resolution gain matrices
%                    will be loaded as well. If false loads only low res.
% OUTPUTS:
%   ConData            - {1 x N_conditions_total} the same cell array as in 
%                        input but some fields are added:
%
%   ConData{c}.G2dLR   - {Nchannels x 2 * NsensorsLR} for low res gain matrices
%                        without dimensionality reduction
%   ConData{c}.G2dHR   - {Nchannels x 2 * NsensorsLR}
%                        the same as above for high res gain matrix
%                        (computed only if bUseHR is set to True)
%   ConData{c}.UP      - {N_sensors_reduced x Nsensors} matrix of left 
%                        singular vectors for G2dLR
%   ConData{c}.G2dLRU  - {Nchannels x N_sensors_reduced} gain matrix 
%                        with reduced number of sensors
% __________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    GainSVDTh = 0.01;
    Nch = length(ChUsed);
    N_conditions_total = length(ConData); % Number of conditions recorded for all subjects altogether
    for c = 1:N_conditions_total
        ConData{c}.NsitesLR = size(ConData{c}.HM_LR.GridLoc, 1);   
        ConData{c}.G2dLR = zeros(Nch, ConData{c}.NsitesLR * 2);
        % reduce tangent space    
        ConData{c}.G2dLR = ReduceTangentSpace(ConData{c}.NsitesLR, ConData{c}.HM_LR.Gain, ChUsed);

        %reduce sensor space
        [ug sg vg] = spm_svd(ConData{c}.G2dLR * ConData{c}.G2dLR', GainSVDTh);
        ConData{c}.UP = ug';
        ConData{c}.G2dLRU = ConData{c}.UP * ConData{c}.G2dLR;
        
        if(bUseHR)
            ConData{c}.NsitesHR = size(ConData{c}.HM_HR.GridLoc, 1);
            ConData{c}.G2dHR = zeros(Nch, ConData{c}.NsitesHR * 2);
            % reduce tangent space
            ConData{c}.G2dHR = ReduceTangentSpace(ConData{c}.NsitesHR, ConData{c}.HM_HR.Gain, ChUsed);
            % Probably need to add here a part for sensor space reduction in case of 
            % high resolution gain matrix
        end;
        c;
    end;