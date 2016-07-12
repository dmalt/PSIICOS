clear all;
setup_subjNames
setup_chLoc

cond_main = '2';
cond_proj_from = '4';
band = [18,22];
tRange = [0,0.3];
sFreq = 500;


% --------------- 3. PSIICOS .3 ------------------ % 
for iSubj = 1:length(subjNames)
% for iSubj = 1:1
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond_main);
	trials_main = LoadTrials(curName, cond_main, band, tRange);
	trials_proj_from = LoadTrials(curName, cond_proj_from, band, tRange);

	CT_m = CrossSpectralTimeseries(trials_main.data);
	CT_m = ProjectAwayFromPowerComplete(CT_m, HM.gain);
	CT_m = RestoreCTdim(CT_m, HM.UP);

	% CT_p = CrossSpectralTimeseries(trials_proj_from.data);
	% CT_p = ProjectAwayFromPowerComplete(CT_p, HM.gain);
	
	% CT_diff = CT_m - CT_p;
	conInds_full{iSubj} = SensorConnectivity((CT_m), 100);
	conInds_real{iSubj} = SensorConnectivity(real(CT_m), 100);
	conInds_imag{iSubj} = SensorConnectivity(imag(CT_m), 100);
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