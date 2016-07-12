clear all;
setup_subjNames;
setup_chLoc;
cond_main = '2';
cond_proj_from = '4';

sFreq = 500;
band_erp = [1,10]; % band for ERP
band_trials = [18,22];
tRange = [-0.5,1];

% ----------------------- Without projection, cond 2  --------------------------- %
for iSubj = 1:length(subjNames)
	%---------- Align trials ------------%
	curName = subjNames{iSubj};
	trials_erp = LoadTrials(curName, cond_main, band_erp, tRange);
	trials_main = LoadTrials(curName, cond_main, band_trials, tRange);

	trials_erp = trials_erp.data;
	trials_data = trials_main.data;
	tr_aligned = alignTrials(trials_erp, trials_data);
	%------------------------------------%

	HM = LoadHeadModel(curName, cond_main);

	trials_proj_from = LoadTrials(curName, cond_proj_from, band, tRange);
	tr_proj_from_aligned = alignTrials(trials_proj_from.data);

	CT_m = CrossSpectralTimeseries(tr_aligned);
	CT_p = CrossSpectralTimeseries(tr_proj_from_aligned);

	CT_m = ProjectAwayFromPowerComplete(CT_m, HM.gain);
	CT_p = ProjectAwayFromPowerComplete(CT_p, HM.gain);

	% CT = ProjFromCond(CT_m, CT_p);
	% CT = RestoreCTdim(CT, HM.UP);


	iRnk = 1;
	for rnk = 0:20
		CT_2_from_4 = ProjFromCond(CT_m, CT_p, rnk);
		CT_2_from_4 = RestoreCTdim(CT_2_from_4, HM.UP);
		conInds_2_from_4{iSubj}{iRnk} = SensorConnectivity((CT_2_from_4), 100);
		iRnk = iRnk + 1;
	end

	% conInds_full{iSubj} = SensorConnectivity((CT), 100);
	% conInds_real{iSubj} = SensorConnectivity(real(CT), 100);
	% conInds_imag{iSubj} = SensorConnectivity(imag(CT), 100);
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
figure;errorbar(0:20, mean_2_from_4, std_2_from_4 ./ sqrt(10), 'bs-')
	% ------------------------------------------------ %


% count = 1;
% for s1 = 1:10
%    for s2 = s1+1:10
%        if(s1~=s2)
%            CS_full(count) = ConnectivitySimilarity(conInds_full{s1}, conInds_full{s2}, ChLoc);
%            CS_real(count) = ConnectivitySimilarity(conInds_real{s1}, conInds_real{s2}, ChLoc);
%            CS_imag(count) = ConnectivitySimilarity(conInds_imag{s1}, conInds_imag{s2}, ChLoc);
%            count = count + 1;
%        end
%    end;
% end;

% disp('Cond 2')
% mean_full = mean(CS_full)
% mean_real = mean(CS_real)
% mean_imag = mean(CS_imag)
% % ------------------------------------------------ %

% --------------------------------------------------------------------- %