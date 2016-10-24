clear conInds_full;
clear conInds_real;
clear conInds_imag;
clear CT_final;

subjNames = {   '0003_pran', ... % 1
				'0019_shev', ... % 2
				'0030_koal', ... % 3
				'0062_peek', ... % 4
				'0074_kuni', ... % 5
				'0100_kase', ... % 6
				'0106_supo', ... % 7
				'0108_bami', ... % 8
				'0109_zvma', ... % 9
				'0130_hagr'};    % 10

sensors_file = '/home/dmalt/ups/data/channel_vectorview306.mat';
ChUsed = ups.PickChannels('grad');
ChLoc = ups.ReadChannelLocations(sensors_file, ChUsed);

main_cond = '4';
proj_cond = '1';
band = [18,22];
tRange = [0,0.7];
sFreq = 500;
d_min = 0.05; % minimal allowed length for connections
nConn_ini = 100; % initial number of connections
nConn = 100; % number of connections to end up with
pwrRnk = 500;
protocol_path = '/home/dmalt/PSIICOS_osadtchii'
% -------------------------- <><><><><><><><><><> ----------------------------- % 
for iSubj = 1:length(subjNames)
	curName = subjNames{iSubj};

	HM = ups.LoadHeadModel(curName, main_cond);

	CT_m = ups.GetCTS(curName, main_cond, band, tRange, 0.01, protocol_path, true, true);
	CT_prf = ups.GetCTS(curName, proj_cond, band, tRange, 0.01, protocol_path, true, true);

	CT_m = ps.ProjectAwayFromPowerComplete(CT_m, HM.gain, pwrRnk);
	CT_prf = ps.ProjectAwayFromPowerComplete(CT_prf, HM.gain, pwrRnk);

	CT_final{iSubj} = ps.ProjFromCond(CT_m, CT_prf, 6);
	% CT_final{iSubj} = CT_m; 

	CT_final{iSubj} = ups.RestoreCTdim(CT_final{iSubj}, HM.UP);
	% CT_final{iSubj} = CT_final{iSubj}(:,201:250);

	conInds_full{iSubj} = ups.GetSensorConnectivity(CT_final{iSubj},       nConn_ini);
	conInds_real{iSubj} = ups.GetSensorConnectivity(real(CT_final{iSubj}), nConn_ini);
	conInds_imag{iSubj} = ups.GetSensorConnectivity(imag(CT_final{iSubj}), nConn_ini);

	% conInds_full{iSubj} = DropLongConn(conInds_full{iSubj}, ChLoc, d_min, nConn);
	% conInds_real{iSubj} = DropLongConn(conInds_real{iSubj}, ChLoc, d_min, nConn);
	% conInds_imag{iSubj} = DropLongConn(conInds_imag{iSubj}, ChLoc, d_min, nConn);
end


conn_metrics_script;
draw_conn_sim;

mean_full
mean_real
mean_imag

p_full_real = ranksum(CS_full, CS_real)
p_full_imag = ranksum(CS_full, CS_imag)
p_imag_real = ranksum(CS_imag, CS_real)
% ------------------------- <><><><><><><><><><><> ---------------------------- %
