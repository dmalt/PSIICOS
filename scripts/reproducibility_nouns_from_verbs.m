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

cond_noun = '2';
cond_verb = '4';
band = [16,25];
tRange = [0,0.7];
sFreq = 500;


% --------------- 3. PSIICOS .3 ------------------ % 
h_fig = figure;
set(h_fig, 'name', 'PSIICOS', 'numbertitle', 'off');
for iSubj = 1:length(subjNames)
% for iSubj = 1:1
	curName = subjNames{iSubj};
	HM = LoadHeadModel(curName, cond_noun);
	trials_noun = LoadTrials(curName, cond_noun, band, tRange);
	trials_verb = LoadTrials(curName, cond_verb, band, tRange);

	CT_noun = CrossSpectralTimeseries(trials_noun.data);
	CT_verb = CrossSpectralTimeseries(trials_verb.data);
	CT = ProjAwayFromCond(CT_noun, CT_verb);

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