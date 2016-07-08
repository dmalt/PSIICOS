clear all;
setup_subjNames
setup_chLoc

cond_noun = '2';
cond_pseudoword = '4';
band = [16,25];
tRange = [0,0.7];
sFreq = 500;


% --------------- 3. PSIICOS .3 ------------------ % 
for iSubj = 1:length(subjNames)
% for iSubj = 1:1
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond_noun);
	trials_2 = LoadTrials(curName, cond_noun, band, tRange);
	trials_4 = LoadTrials(curName, cond_pseudoword, band, tRange);

	CT_2 = CrossSpectralTimeseries(trials_2.data);
	CT_2 = ProjectAwayFromPowerComplete(CT_2, HM.gain);

	CT_4 = CrossSpectralTimeseries(trials_4.data);
	CT_4 = ProjectAwayFromPowerComplete(CT_4, HM.gain);
	iRnk = 1;
	
	CT_diff = CT_2 - CT_4;
	CT_diff = RestoreCTdim(CT_diff, HM.UP);
	conInds_diff{iSubj} = SensorConnectivity((CT_diff), 100);
end


count = 1;
for s1 = 1:10
   for s2 = s1+1:10
       if(s1~=s2)
           CS_diff(count) = ConnectivitySimilarity(conInds_diff{s1}, conInds_diff{s2}, ChLoc);
           count = count + 1;
       end
   end;
end;


mean_diff = mean(CS_diff)

	% ------------------------------------------------ %