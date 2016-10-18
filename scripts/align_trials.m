clear all;
setup_subjNames;
setup_chLoc;
cond_main = '4';
cond_proj_from = '4';

sFreq = 500;
band_erp = [1,10]; % band for ERP
band_tr = [18,22]; % band for ERP
tRange = [-0.5,1];

% ----------------------- Without projection, cond 2  --------------------------- %
for iSubj = 1:length(subjNames)
	%---------- Align trials ------------%
	curName = subjNames{iSubj};
	trials_erp = LoadTrials(curName, cond_main, band_erp, tRange);
	trials_main = LoadTrials(curName, cond_main, band_tr, tRange);

	trials_data = trials_main.data;
	trials_data_erp = trials_erp.data;

	tr_aligned = alignTrials(trials_data_erp, trials_data);
	%------------------------------------%

	HM = LoadHeadModel(curName, cond_main);

	CT_m = CrossSpectralTimeseries(tr_aligned);
	CT_m = ProjectAwayFromPowerComplete(CT_m, HM.gain);
	CT_m = RestoreCTdim(CT_m, HM.UP);

	conInds_full{iSubj} = GetSensorConnectivity((CT_m), 100);
	conInds_real{iSubj} = GetSensorConnectivity(real(CT_m), 100);
	conInds_imag{iSubj} = GetSensorConnectivity(imag(CT_m), 100);
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

disp('Cond 2')
mean_full = mean(CS_full)
mean_real = mean(CS_real)
mean_imag = mean(CS_imag)
% ------------------------------------------------ %

% --------------------------------------------------------------------- %