clear all;
setup_subjNames
setup_chLoc

cond_noun = '4';
cond_pseudoword = '1';
band = [18,22];
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

	CT_2_from_4 = ProjFromCond(CT_2, CT_4, 6);
	CT_2_from_4 = RestoreCTdim(CT_2_from_4, HM.UP);
	% conInds_2_from_4{iSubj} = SensorConnectivity((CT_2_from_4), 100);

	conInds_full{iSubj} = SensorConnectivity((CT_2_from_4), 100);
    conInds_real{iSubj} = SensorConnectivity(real(CT_2_from_4), 100);
    conInds_imag{iSubj} = SensorConnectivity(imag(CT_2_from_4), 100);

end

count = 1;
for s1 = 1:10
   for s2 = s1+1:10
       if(s1~=s2)
           CS_full(count) = ConnectivitySimilarity(conInds_full{s1}, conInds_full{s2}, ChLoc);
           CS_real(count) = ConnectivitySimilarity(conInds_real{s1}, conInds_real{s2}, ChLoc);
           CS_imag(count) = ConnectivitySimilarity(conInds_imag{s1}, conInds_imag{s2}, ChLoc);
           count = count + 1;
       end
   end;
end;


mean_full = mean(CS_full)
mean_real = mean(CS_real)
mean_imag = mean(CS_imag)
	% ------------------------------------------------ %