clear all;
setup_subjNames
setup_chLoc

cond_noun = '2';
cond_verb = '4';
band = [16,25];
tRange = [0,0.7];
sFreq = 500;


% --------------- 3. PSIICOS .3 ------------------ % 
for iSubj = 1:length(subjNames)
% for iSubj = 1:1
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond_noun);
	trials_2 = LoadTrials(curName, cond_noun, band, tRange);
	trials_4 = LoadTrials(curName, cond_verb, band, tRange);

	CT_2 = CrossSpectralTimeseries(trials_2.data);
	CT_2 = ProjectAwayFromPowerComplete(CT_2, HM.gain);

	CT_4 = CrossSpectralTimeseries(trials_4.data);
	CT_4 = ProjectAwayFromPowerComplete(CT_4, HM.gain);
	iRnk = 1;
	for rnk = 0:20
		CT_2_from_4 = ProjFromCond(CT_2, CT_4, rnk);
		CT_2_from_4 = RestoreCTdim(CT_2_from_4, HM.UP);
		conInds_2_from_4{iSubj}{iRnk} = SensorConnectivity((CT_2_from_4), 100);
		iRnk = iRnk + 1;
	end
end

iRnk = 1;
for rnk = 0:20
	count = 1;
	for s1 = 1:10
	   for s2 = s1+1:10
	       if(s1~=s2)
	           CS_2_from_4(count, iRnk) = ConnectivitySimilarity(conInds_2_from_4{s1}{iRnk}, conInds_2_from_4{s2}{iRnk}, ChLoc);
	           count = count + 1;
	       end
	   end;
	end;

	mean_2_from_4(iRnk) = mean(CS_2_from_4(:, iRnk));
	std_2_from_4(iRnk) = std(CS_2_from_4(:,iRnk));
	iRnk = iRnk + 1;
end
errorbar(0:20, mean_2_from_4, std_2_from_4, 'bs-')
	% ------------------------------------------------ %