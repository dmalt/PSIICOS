% clear all;
setup_subjNames
setup_chLoc

main_cond = '2';
proj_cond = '1';
band = [18,22];
tRange = [0,0.7];
sFreq = 500;
d_min = 0.05; % minimal allowed length for connections
nConn_ini = 200; % initial number of connections
nConn = 100; % number of connections to end up with

% --------------- 3. PSIICOS .3 ------------------ % 
for iSubj = 1:length(subjNames)
	curName = subjNames{iSubj};

	HM = LoadHeadModel(curName, main_cond);

	tr_m   = LoadTrials(curName, main_cond, band, tRange);
	tr_prf = LoadTrials(curName, proj_cond, band, tRange);

	CT_m   = CrossSpectralTimeseries(tr_m.data);
	CT_prf = CrossSpectralTimeseries(tr_prf.data);

	CT_m   = ProjectAwayFromPowerComplete(CT_m, HM.gain);
	CT_prf = ProjectAwayFromPowerComplete(CT_prf, HM.gain);

	CT_pr = ProjFromCond(CT_m, CT_prf, 6);

	CT_pr = RestoreCTdim(CT_pr, HM.UP);

	conInds_full{iSubj} = GetSensorConnectivity(CT_pr,       nConn_ini);
	conInds_real{iSubj} = GetSensorConnectivity(real(CT_pr), nConn_ini);
	conInds_imag{iSubj} = GetSensorConnectivity(imag(CT_pr), nConn_ini);

	conInds_full{iSubj} = DropLongConn(conInds_full{iSubj}, ChLoc, d_min, nConn);
	conInds_real{iSubj} = DropLongConn(conInds_real{iSubj}, ChLoc, d_min, nConn);
	conInds_imag{iSubj} = DropLongConn(conInds_imag{iSubj}, ChLoc, d_min, nConn);
end

[mean_full, std_full, CS_full] = ConnSimMetrics(conInds_full, ChLoc);
[mean_real, std_real, CS_real] = ConnSimMetrics(conInds_real, ChLoc);
[mean_imag, std_imag, CS_imag] = ConnSimMetrics(conInds_imag, ChLoc);
	% ------------------------------------------------ %