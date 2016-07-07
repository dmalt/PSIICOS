subjNames = {   '0003_pran', ... % 1
				'0019_shev', ... % 2
				'0030_koal', ... % 3
				'0108_bami', ... % 3
				'0062_peek', ... % 5
				'0074_kuni', ... % 6
				'0100_kase', ... % 7
				'0106_supo', ... % 8
				'0109_zvma', ... % 9
				'0130_hagr'};    % 10
cond = '2';
band = [16,25];
tRange = [0,0.7];
sFreq = 500;
ChUsed = PickChannels('grad');
sensors_file = '~/ups/channel_vectorview306.mat';
ChLoc = ReadChannelLocations(sensors_file, ChUsed);
% ------------- 1. Cross-spectrum .1 ------------ % 
h_fig = figure;
set(h_fig, 'name', 'Cross-spectrum', 'numbertitle', 'off');
for iSubj = 1:length(subjNames)
% for iSubj = 1:1
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond);
	trials = LoadTrials(curName, cond, band, tRange);
	CT = CrossSpectralTimeseries(trials.data);
	CT = RestoreCTdim(CT, HM.UP);
	conInds = SensorConnectivity((CT), 100);
	subplot(2,5,iSubj);
	h_ax = DrawConnectionsOnSensors(conInds);
	% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
	view(h_ax, [0,90]);
end
% ------------------------------------------------ %


% ------- 2. Imaginary cross-spectrum .2 --------- % 
h_fig = figure;
set(h_fig, 'name', 'Imaginary cross-spectrum', 'numbertitle', 'off');
for iSubj = 1:length(subjNames)
% for iSubj = 1:1
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond);
	trials = LoadTrials(curName, cond, band, tRange);
	CT = CrossSpectralTimeseries(trials.data);
	CT = RestoreCTdim(CT, HM.UP);
	conInds = SensorConnectivity(imag(CT), 100);
	subplot(2,5,iSubj);
	h_ax = DrawConnectionsOnSensors(conInds);
	% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
	view(h_ax, [0,90]);
end
% ------------------------------------------------ %


% --------------- 3. PSIICOS .3 ------------------ % 
h_fig = figure;
set(h_fig, 'name', 'PSIICOS', 'numbertitle', 'off');
conInds = {}
for iSubj = 1:length(subjNames)
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond);
	trials = LoadTrials(curName, cond, band, tRange);
	CT = CrossSpectralTimeseries(trials.data);
	CT = ProjectAwayFromPowerComplete(CT, HM.gain);
	CT = RestoreCTdim(CT, HM.UP);
	conInds{iSubj} = SensorConnectivity((CT), 100);
	subplot(2,5,iSubj);
	h_ax = DrawConnectionsOnSensors(conInds{iSubj});
	% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
	view(h_ax, [0,90]);
end


CSa = zeros(10,10);
for s1 = 1:10
   for s2 = 1:10
       if(s1~=s2)
           CSa(s1,s2) = ConnectivitySimilarity(conInds{s1}, conInds{s2}, ChLoc);
       end
   end;
end;
% ------------------------------------------------ %

% ------------- 4. Real PSIICOS .4  -------------- % 
h_fig = figure;
set(h_fig, 'name', 'Real PSIICOS', 'numbertitle', 'off');
for iSubj = 1:length(subjNames)
% for iSubj = 1:1
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond);
	trials = LoadTrials(curName, cond, band, tRange);
	CT = CrossSpectralTimeseries(trials.data);
	CT = ProjectAwayFromPowerComplete(CT, HM.gain);
	CT = RestoreCTdim(CT, HM.UP);
	conInds = SensorConnectivity(real(CT), 100);
	subplot(2,5,iSubj);
	h_ax = DrawConnectionsOnSensors(conInds);
	% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
	view(h_ax, [0,90]);
end
% ------------------------------------------------ %


% ----------- 5. Imaginary PSIICOS .5 ------------ % 
h_fig = figure;
set(h_fig, 'name', 'Imaginary PSIICOS', 'numbertitle', 'off');
for iSubj = 1:length(subjNames)
% for iSubj = 1:1
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond);
	trials = LoadTrials(curName, cond, band, tRange);
	CT = CrossSpectralTimeseries(trials.data);
	CT = ProjectAwayFromPowerComplete(CT, HM.gain);
	CT = RestoreCTdim(CT, HM.UP);
	conInds = SensorConnectivity(imag(CT), 100);
	subplot(2,5,iSubj);
	h_ax = DrawConnectionsOnSensors(conInds);
	% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
	view(h_ax, [0,90]);
end
% ------------------------------------------------ %


% --------------- 6. Coherence .6 ---------------- % 
h_fig = figure;
set(h_fig, 'name', 'Coherence', 'numbertitle', 'off');
for iSubj = 1:length(subjNames)
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond);
	trials = LoadTrials(curName, cond, band, tRange);
	CT = CrossSpectralTimeseries(trials.data);
	CT = RestoreCTdim(CT, HM.UP);
	CohTS = Cp2Coh(CT);
	conInds = SensorConnectivity(CohTS, 100);
	subplot(2,5,iSubj);
	h_ax = DrawConnectionsOnSensors(conInds);
	% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
	view(h_ax, [0,90]);
end
% ------------------------------------------------ %


% ---------- 7. Imaginary coherence .7 ----------- % 
h_fig = figure;
set(h_fig, 'name', 'Imaginary coherence', 'numbertitle', 'off');
for iSubj = 1:length(subjNames)
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond);
	trials = LoadTrials(curName, cond, band, tRange);
	CT = CrossSpectralTimeseries(trials.data);
	CT = RestoreCTdim(CT, HM.UP);
	CohTS = Cp2Coh(CT);
	conInds = SensorConnectivity(imag(CohTS), 100);
	subplot(2,5,iSubj);
	h_ax = DrawConnectionsOnSensors(conInds);
	% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
	view(h_ax, [0,90]);
end
% ------------------------------------------------ %


% ------------------ 8. PLV .8 ------------------- % 
h_fig = figure;
set(h_fig, 'name', 'PLV', 'numbertitle', 'off');
for iSubj = 1:length(subjNames)
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond);
	trials = LoadTrials(curName, cond, band, tRange);
	Tr = RestoreTrDim(trials.data, HM.UP);
	PLV = PLVMatrix(Tr, band, sFreq, false);
	conInds = SensorConnectivity(PLV(:), 100);
	subplot(2,5,iSubj);
	h_ax = DrawConnectionsOnSensors(conInds);
	% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
	view(h_ax, [0,90]);
end
% ------------------------------------------------ %


% ------------------ 9. PLI .9 ------------------- % 
h_fig = figure;
set(h_fig, 'name', 'PLI', 'numbertitle', 'off');
for iSubj = 1:length(subjNames)
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond);
	trials = LoadTrials(curName, cond, band, tRange);
	Tr = RestoreTrDim(trials.data, HM.UP);
	PLI = PLIMatrix(Tr, band, sFreq, false);
	conInds = SensorConnectivity(PLI(:), 100);
	subplot(2,5,iSubj);
	h_ax = DrawConnectionsOnSensors(conInds);
	% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
	view(h_ax, [0,90]);
end
% ------------------------------------------------ %

% ------------------ 10. wPLI .10 ------------------ % 
h_fig = figure;
set(h_fig, 'name', 'wPLI', 'numbertitle', 'off');
for iSubj = 1:length(subjNames)
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond);
	trials = LoadTrials(curName, cond, band, tRange);
	Tr = RestoreTrDim(trials.data, HM.UP);
	wPLI = wPLIMatrix(Tr, band, sFreq, false);
	conInds = SensorConnectivity(wPLI(:), 100);
	subplot(2,5,iSubj);
	h_ax = DrawConnectionsOnSensors(conInds);
	% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
	view(h_ax, [0,90]);
end
% ------------------------------------------------ %

