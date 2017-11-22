clear all;


protocol_path = '/home/dmalt/fif_matlab/Brainstorm_db/motorprobe';
subj_ID = 'AP'
suffix = '_control_RH_raw_tsss_ica';
condition = [subj_ID, suffix];
freq_band = [10, 20];
tRange = [-0.9, -0.5];
isLR = true;
GainSVDTh = 0.01;
sFreq = 500;

ChUsed = ups.PickChannels('grad');
sensors_file = '~/ups/data/channel_vectorview306.mat';
ChLoc = ups.ReadChannelLocations(sensors_file, ChUsed);
% ------------- 1. Cross-spectrum .1 ------------ % 
h_fig = figure;
set(h_fig, 'name', 'Cross-spectrum', 'numbertitle', 'off');

HM = ups.LoadHeadModel(subj_ID, condition, protocol_path, isLR, GainSVDTh);
trials = ups.LoadTrials(subj_ID, condition, freq_band, tRange, GainSVDTh, protocol_path);
CT = ups.CrossSpectralTimeseries(trials.data);
CT = ups.RestoreCTdim(CT, HM.UP);
conInds = ups.GetSensorConnectivity(real(CT), 100);
h_ax = ups.DrawConnectionsOnSensors(conInds);
% h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
view(h_ax, [0,90]);
% ------------------------------------------------ %


% % ------- 2. Imaginary cross-spectrum .2 --------- % 
% h_fig = figure;
% set(h_fig, 'name', 'Imaginary cross-spectrum', 'numbertitle', 'off');
% for iSubj = 1:length(subjNames)imag
% % for iSubj = 1:1
%   subj_ID = subjNames{iSubj};
%   HM = LoadHeadModel(subj_ID, cond);
%   trials = LoadTrials(subj_ID, cond, band, tRange);
%   CT = ups.CrossSpectralTimeseries(trials.data);
%   CT = RestoreCTdim(CT, HM.UP);
%   conInds = GetSensorConnectivity(imag(CT), 100);
%   subplot(2,5,iSubj);
%   h_ax = DrawConnectionsOnSensors(conInds);
%   % h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
%   view(h_ax, [0,90]);
% end
% % ------------------------------------------------ %


% % --------------- 3. PSIICOS .3 ------------------ % 
% h_fig = figure;
% set(h_fig, 'name', 'PSIICOS', 'numbertitle', 'off');
% conInds = {}
% for iSubj = 1:length(subjNames)
%   subj_ID = subjNames{iSubj};
%   HM = LoadHeadModel(subj_ID, cond);
%   trials = LoadTrials(subj_ID, cond, band, tRange);
%   CT = ups.CrossSpectralTimeseries(trials.data);
%   CT = ProjectAwayFromPowerComplete(CT, HM.gain);
%   CT = RestoreCTdim(CT, HM.UP);
%   conInds{iSubj} = GetSensorConnectivity((CT), 100);
%   subplot(2,5,iSubj);
%   h_ax = DrawConnectionsOnSensors(conInds{iSubj});
%   % h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
%   view(h_ax, [0,90]);
% end


% CSa = zeros(10,10);
% for s1 = 1:10
%    for s2 = 1:10
%        if(s1~=s2)
%            CSa(s1,s2) = ConnectivitySimilarity(conInds{s1}, conInds{s2}, ChLoc);
%        end
%    end;
% end;
% % ------------------------------------------------ %

% % ------------- 4. Real PSIICOS .4  -------------- % 
% h_fig = figure;
% set(h_fig, 'name', 'Real PSIICOS', 'numbertitle', 'off');
% for iSubj = 1:length(subjNames)
% % for iSubj = 1:1
%   subj_ID = subjNames{iSubj};
%   HM = LoadHeadModel(subj_ID, cond);
%   trials = LoadTrials(subj_ID, cond, band, tRange);
%   CT = ups.CrossSpectralTimeseries(trials.data);
%   CT = ProjectAwayFromPowerComplete(CT, HM.gain);
%   CT = RestoreCTdim(CT, HM.UP);
%   conInds = GetSensorConnectivity(real(CT), 100);
%   subplot(2,5,iSubj);
%   h_ax = DrawConnectionsOnSensors(conInds);
%   % h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
%   view(h_ax, [0,90]);
% end
% % ------------------------------------------------ %


% % ----------- 5. Imaginary PSIICOS .5 ------------ % 
% h_fig = figure;
% set(h_fig, 'name', 'Imaginary PSIICOS', 'numbertitle', 'off');
% for iSubj = 1:length(subjNames)
% % for iSubj = 1:1
%   subj_ID = subjNames{iSubj};
%   HM = LoadHeadModel(subj_ID, cond);
%   trials = LoadTrials(subj_ID, cond, band, tRange);
%   CT = ups.CrossSpectralTimeseries(trials.data);
%   CT = ProjectAwayFromPowerComplete(CT, HM.gain);
%   CT = RestoreCTdim(CT, HM.UP);
%   conInds = GetSensorConnectivity(imag(CT), 100);
%   subplot(2,5,iSubj);
%   h_ax = DrawConnectionsOnSensors(conInds);
%   % h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
%   view(h_ax, [0,90]);
% end
% % ------------------------------------------------ %


% % --------------- 6. Coherence .6 ---------------- % 
% h_fig = figure;
% set(h_fig, 'name', 'Coherence', 'numbertitle', 'off');
% for iSubj = 1:length(subjNames)
%   subj_ID = subjNames{iSubj};
%   HM = LoadHeadModel(subj_ID, cond);
%   trials = LoadTrials(subj_ID, cond, band, tRange);
%   CT = ups.CrossSpectralTimeseries(trials.data);
%   CT = RestoreCTdim(CT, HM.UP);
%   CohTS = Cp2Coh(CT);
%   conInds = GetSensorConnectivity(CohTS, 100);
%   subplot(2,5,iSubj);
%   h_ax = DrawConnectionsOnSensors(conInds);
%   % h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
%   view(h_ax, [0,90]);
% end
% % ------------------------------------------------ %


% % ---------- 7. Imaginary coherence .7 ----------- % 
% h_fig = figure;
% set(h_fig, 'name', 'Imaginary coherence', 'numbertitle', 'off');
% for iSubj = 1:length(subjNames)
%   subj_ID = subjNames{iSubj};
%   HM = LoadHeadModel(subj_ID, cond);
%   trials = LoadTrials(subj_ID, cond, band, tRange);
%   CT = ups.CrossSpectralTimeseries(trials.data);
%   CT = RestoreCTdim(CT, HM.UP);
%   CohTS = Cp2Coh(CT);
%   conInds = GetSensorConnectivity(imag(CohTS), 100);
%   subplot(2,5,iSubj);
%   h_ax = DrawConnectionsOnSensors(conInds);
%   % h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
%   view(h_ax, [0,90]);
% end
% % ------------------------------------------------ %


% % ------------------ 8. PLV .8 ------------------- % 
% h_fig = figure;
% set(h_fig, 'name', 'PLV', 'numbertitle', 'off');
% for iSubj = 1:length(subjNames)
%   subj_ID = subjNames{iSubj};
%   HM = LoadHeadModel(subj_ID, cond);
%   trials = LoadTrials(subj_ID, cond, band, tRange);
%   Tr = RestoreTrDim(trials.data, HM.UP);
%   PLV = PLVMatrix(Tr, band, sFreq, false);
%   conInds = GetSensorConnectivity(PLV(:), 100);
%   subplot(2,5,iSubj);
%   h_ax = DrawConnectionsOnSensors(conInds);
%   % h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
%   view(h_ax, [0,90]);
% end
% % ------------------------------------------------ %


% % ------------------ 9. PLI .9 ------------------- % 
% h_fig = figure;
% set(h_fig, 'name', 'PLI', 'numbertitle', 'off');
% for iSubj = 1:length(subjNames)
%   subj_ID = subjNames{iSubj};
%   HM = LoadHeadModel(subj_ID, cond);
%   trials = LoadTrials(subj_ID, cond, band, tRange);
%   Tr = RestoreTrDim(trials.data, HM.UP);
%   PLI = PLIMatrix(Tr, band, sFreq, false);
%   conInds = GetSensorConnectivity(PLI(:), 100);
%   subplot(2,5,iSubj);
%   h_ax = DrawConnectionsOnSensors(conInds);
%   % h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
%   view(h_ax, [0,90]);
% end
% % ------------------------------------------------ %

% % ------------------ 10. wPLI .10 ------------------ % 
% h_fig = figure;
% set(h_fig, 'name', 'wPLI', 'numbertitle', 'off');
% for iSubj = 1:length(subjNames)
%   subj_ID = subjNames{iSubj};
%   HM = LoadHeadModel(subj_ID, cond);
%   trials = LoadTrials(subj_ID, cond, band, tRange);
%   Tr = RestoreTrDim(trials.data, HM.UP);
%   wPLI = wPLIMatrix(Tr, band, sFreq, false);
%   conInds = GetSensorConnectivity(wPLI(:), 100);
%   subplot(2,5,iSubj);
%   h_ax = DrawConnectionsOnSensors(conInds);
%   % h_title = title(h_ax, [num2str(band), ' Hz ', num2str(tRange), ' s'] );
%   view(h_ax, [0,90]);
% end
% % ------------------------------------------------ %

