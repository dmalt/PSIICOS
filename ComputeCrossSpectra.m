function ConData = ComputeCrossSpectra(ConData)
    for sc = 1:length(ConData)
        fprintf('%d Computing cross-spectral matrix ....\n' , sc); 
        ConData{sc}.CrossSpecTime = CrossSpectralTimeseries(ConData{sc}.Trials); 
        ConData{sc}.CrossSpecTimeInd = CrossSpectralTimeseries(ConData{sc}.Trials,true);
        % compute their projected versions                                                      
        [ConData{sc}.CrossSpecTimeP, ConData{sc}.Upwr] = ProjectAwayFromPowerFixedOr(ConData{sc}.CrossSpecTime, ConData{sc}.G2dLRU,350);

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
