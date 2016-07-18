clear conInds_full;
clear conInds_real;
clear conInds_imag;
setup_subjNames
setup_chLoc

main_cond = '2';
proj_cond = '1';
band = [18,22];
tRange = [0.2,0.3];
sFreq = 500;
d_min = 0.05; % minimal allowed length for connections
nConn_ini = 200; % initial number of connections
nConn = 100; % number of connections to end up with

% -------------------------- <><><><><><><><><><> ----------------------------- % 
for iSubj = 1:length(subjNames)
	curName = subjNames{iSubj};

	HM = LoadHeadModel(curName, main_cond);

	% tr_m   = LoadTrials(curName, main_cond, band, tRange);
	% tr_prf = LoadTrials(curName, proj_cond, band, tRange);

	CT_m = GetCTS(curName, main_cond, band, tRange);
	CT_prf = GetCTS(curName, proj_cond, band, tRange);

	% CT_m   = CrossSpectralTimeseries(tr_m.data);
	% CT_prf = CrossSpectralTimeseries(tr_prf.data);

	CT_m   = ProjectAwayFromPowerComplete(CT_m, HM.gain);
	CT_prf = ProjectAwayFromPowerComplete(CT_prf, HM.gain);

	CT_final = ProjFromCond(CT_m, CT_prf, 6);
	% CT_final = CT_m; 

	CT_final = RestoreCTdim(CT_final, HM.UP);

	conInds_full{iSubj} = GetSensorConnectivity(CT_final,       nConn_ini);
	conInds_real{iSubj} = GetSensorConnectivity(real(CT_final), nConn_ini);
	conInds_imag{iSubj} = GetSensorConnectivity(imag(CT_final), nConn_ini);

	conInds_full{iSubj} = DropLongConn(conInds_full{iSubj}, ChLoc, d_min, nConn);
	conInds_real{iSubj} = DropLongConn(conInds_real{iSubj}, ChLoc, d_min, nConn);
	conInds_imag{iSubj} = DropLongConn(conInds_imag{iSubj}, ChLoc, d_min, nConn);
end

[mean_full, std_full, CS_full] = ConnSimMetrics(conInds_full, ChLoc);
[mean_real, std_real, CS_real] = ConnSimMetrics(conInds_real, ChLoc);
[mean_imag, std_imag, CS_imag] = ConnSimMetrics(conInds_imag, ChLoc);

mean_full
mean_real
mean_imag

p_full_real = ranksum(CS_full, CS_real)
p_full_imag = ranksum(CS_full, CS_imag)
p_imag_real = ranksum(CS_imag, CS_real)
% ------------------------- <><><><><><><><><><><> ---------------------------- %