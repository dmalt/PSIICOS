[HM, CT, Tr, Ctx] = GenerData(pi/2, 100); %, 0.01, 0.35, 0, true);
% CTp = ProjectAwayFromPowerComplete(CT, HM.gain);
CT = RestoreCTdim(CT, HM.UP);

% CohTS = Cp2Coh(CT);
% ConInds = GetSensorConnectivity(imag(CohTS), 20);
% ax = DrawConnectionsOnSensors(ConInds);
% title(ax, 'ImCoh');
% nTr = size(Tr,3);

% for iTr = 1:nTr
% 	TrRestored(:,:,iTr) = HM.UP' * squeeze(Tr(:,:,iTr));
% end

TrRestored = RestoreTrDim(Tr, HM.UP);

PLV = wPLIMatrix(TrRestored, [9.5,10.5], 500);
conInds = GetSensorConnectivity(PLV(:), 100);
DrawConnectionsOnSensors(conInds);
% CTp = RestoreCTdim(CTp, HM.UP);
% conInds = GetSensorConnectivity(imag(CT), 0.3);
% DrawConnectionsOnSensors(conInds);
% conIndsp = GetSensorConnectivity(CTp, 0.3);
% DrawConnectionsOnSensors(conIndsp);