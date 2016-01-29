function ConData = ComputeCrossSpectra(ConData)
% ------------------------------------------------------------------
% ComputeCrossSpectra: compute cross-spectral matrices based
% on ConData{c}.Trials with and without projection from volume
% conduction operation for induced activity and per se 
% and store the result in the same cell array
% ------------------------------------------------------------------
% FORMAT:
%   ConData = ComputeCrossSpectra(ConData) 
% INPUTS:
%   ConData        - {1 x N_conditions_total} cell array with each
%                    element containing a structure for a single
%                    condition recording from protocol.
% OUTPUTS:
%   ConData:
%   ConData{c}.CrossSpecTime       - {N_sensors_reduced ^ 2 x Ntimes}
%                                    cross-spectum matrix for sensors
%   ConData{c}.CrossSpecTimeInd    - {N_sensors_reduced ^ 2 x Ntimes}
%                                    cross-spectrum matrix for induced
%                                    activity on sensors
%   ConData{c}.CrossSpecTimeP      - {N_sensors_reduced ^ 2 x Ntimes}                               
%                                    projected from VC cross-spectrum 
%                                    matrix on sensors 
%   ConData{c}.CrossSpecTimeIndP   - {N_sensors_reduced ^ 2 x Ntimes}
%                                    projected from VC cross-spectrum
%                                    matrix for induced activity on 
%                                    sensors
% ________________________________________________________________________
% Alex Ossadtchii ossadtchi@gmail.com, Dmitrii Altukhov, dm.altukhov@ya.ru

    for sc = 1:length(ConData)
        fprintf('%d Computing cross-spectral matrix ....\n' , sc); 
        ConData{sc}.CrossSpecTime = CrossSpectralTimeseries(ConData{sc}.Trials); 
        ConData{sc}.CrossSpecTimeInd = CrossSpectralTimeseries(ConData{sc}.Trials, true);
        % compute their projected versions                                                      
        [ConData{sc}.CrossSpecTimeP, ConData{sc}.Upwr] = ...
                ProjectAwayFromPowerComplete(ConData{sc}.CrossSpecTime, ConData{sc}.G2dLRU, 350);

        ConData{sc}.CrossSpecTimeIndP = ConData{sc}.CrossSpecTimeInd - ...
                                            ConData{sc}.Upwr * ConData{sc}.Upwr' * ...
                                            ConData{sc}.CrossSpecTimeInd;
        %UP
        % if(bComputePLI)
        %     Trials = zeros(size(ConData{sc}.UP, 2), size(ConData{sc}.Trials, 2), size(ConData{sc}.Trials, 2));
        %     for tr = 1:size(ConData{sc}.Trials, 3)
        %         Trials(:,:,tr) = ConData{sc}.UP' * ConData{sc}.Trials(:,:,tr);
        %     end;
        %     ConData{sc}.wPLI =  wPLIMatrix(Trials(:,1:256,:), Band, Fsamp, true);
        % end;
    end;
