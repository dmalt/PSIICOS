%clear all;
close all;
%% Parameters block
InducedScale = 0.2;
EvokedScale = 0.0/2.5;
GainSVDTh = 0.001; % 0.05 results into 47 eigensensors and makes it run faster but produces less contrasting subcorr scans
% for a more reliable preformance use 0.01 to get all the sensor on board but be ready to wait;
NetworkPairIndex{1} = [1,2];
NetworkPairIndex{2} = [1,2,3];
%% Load forward model and reduce it  
% load reduced forward model (GLowRes)
load('GLowRes.mat'); 
% get grid node locations
R = GLowRes.GridLoc;
% set to use gradiometers only
ChUsed = 1:306; 
ChUsed(3:3:end) = [];
% calculate tangential plane dipoles
[Nch, Nsites] = size(GLowRes.Gain(ChUsed,1:3:end));
G2d = zeros(Nch,Nsites*2);
G2d0 = zeros(Nch,Nsites*2);
range = 1:2;
for i=1:Nsites
    g = [GLowRes.Gain(ChUsed,1+3*(i-1)) GLowRes.Gain(ChUsed,2+3*(i-1)) GLowRes.Gain(ChUsed,3+3*(i-1))];
    [u sv v] = svd(g);
    gt = g*v(:,1:2);
      G2d(:,range) = gt*diag(1./sqrt(sum(gt.^2,1)));
      G2d0(:,range) = gt;
    range = range + 2;
end;
% reduce the sensor space
[ug sg vg] = spm_svd(G2d*G2d',GainSVDTh);
UP = ug';
G2dU = UP*G2d;
G2d0U = UP*G2d0;

band = [8 12];
PHI = [pi/20 pi/2-pi/20];
it = 1;
for phi = PHI 
    %% Data simulation
    if(it==1)
        % we will use the same phase shifts in the second and subsequent
        % iterations. We will store the random phases in PhaseShifts 
        [ Evoked, Induced, BrainNoise, SensorNoise, G2dHiRes, RHiRes, Fs,Ntr, XYZGen, Ggen, PhaseShifts] = ...
        SimulateDataPhase(NetworkPairIndex{2},phi,false,[]);
    else
        % and use PhaseShits from the first iteration
        [ Evoked, Induced, BrainNoise, SensorNoise, G2dHiRes, RHiRes, Fs,Ntr, XYZGen, Ggen, PhaseShifts] = ...
        SimulateDataPhase(NetworkPairIndex{2},phi,false,PhaseShifts);
    end;        
    % mix noise and data 
    % in order to control SNR we first normalize norm(BrainNoise(:)) = 1 and 
    % norm(Induced(:)) = 1 and then mix the two with the coefficient
    Data0 = BrainNoise + InducedScale*Induced + EvokedScale*Evoked ;
   % Data0 = InducedScale*Induced + EvokedScale*Evoked ;
    [bf af] = butter(5,band/(0.5*Fs),'bandpass');
    % Filter in the band of interest
    Data = filtfilt(bf,af,Data0')';
    %% Reshape the data in a 3D structure(Space x Time x Epochs)
    [Nch Tcnt] = size(Data);
    T = fix(Tcnt/Ntr);
    Nch = size(UP,1);
    %reshape Data and store in a 3D array X
    X1 = zeros(Nch,T,Ntr);
    X0 = zeros(size(Data,1),T,Ntr);
    range = 1:T;
    for i=1:Ntr
        X1(:,:,i) = UP*Data(:,range);
        X0(:,:,i) = Data(:,range);
        range = range+T;
    end;
    %calculate SNR
    NoiseF  = filtfilt(bf,af,(UP*BrainNoise)')';
    SignalF  = filtfilt(bf,af,(InducedScale*UP*Induced)')';
    SNR{it} = norm(SignalF(:))/norm(NoiseF(:));
    
    %% Calculate band cross-spectral matrix 
    CrossSpecTime{it} = CrossSpectralTimeseries(X1);
    C{it} = reshape(mean(CrossSpecTime{it},2),Nch,Nch);
   % [Cfft{it} Cohfft{it}] = CrossSpectralMatrix(X1,[9.5 10.5],Fs,false);
   % [Cfftvec] = ProjectAwayFromPowerFixedOr(Cfft{it}(:), G2dU);
   % Cfftp{it} = reshape(Cfftvec,size(Cfft{it}));
    % because of nonlinearity use original sensor space data to compute the
    % PLV
   % [PLV{it} ] = PLVMatrix(X0,[9.5 10.5],Fs,false);
   % [PLI{it} ] = PLIMatrix(X0,[9.5 10.5],Fs,false);

    
    %% Apply PSIICOS
   % [Qpsiicos{it}, IND{it}, Cp{it}, Upwr{it}] = PSIICOS_Fast(C{it}, G2dU,350);
    [ indep_topo, c_ss_hat, PVU, SCOR, INDrap, Cp, Upwr]=RAP_PSIICOS_Fast(C{it}, G2dU,3);
    [SPCpsiicos{it}, TPRpsiicos{it}] = GenerateROC(Qpsiicos{it},0.015,R, IND{it}, 100,XYZGen,NetworkPairIndex{2});
    
    %% Apply iDICS
    [Qdics{it}, Psdics{it}, IND{it}] = iDICS(C{it},G2dU);
    [SPCdics{it}, TPRdics{it}] = GenerateROC(Qdics{it},0.015,R, IND{it}, 100,XYZGen,NetworkPairIndex{2});
    it = it+1;
end;
%% Make a figure
Labels = {'r.-','b.-','ro-','bo-','m.'};
figure
for it = 1:length(PHI)
    plot(1-SPCpsiicos{it},TPRpsiicos{it},Labels{it*length(PHI)-1});
    hold on;
    plot(1-SPCdics{it},TPRdics{it},Labels{it*length(PHI)});
end
axis([0 0.01 0 1.1])
grid
legend('PSIICOS, phi = pi/10','iDICS, phi = pi/10', 'PSIICOS, phi = pi/2','iDICS, phi = pi/2');
xlabel('1-specificity');
ylabel('sensitivity');

% 
% load('c:\MyWriteups\iRAPMusicPaper\Simulations\MEGSensors.mat');
% for i = 1:length(ChUsed)
%     ChLoc(:,i) = MEGSensors.Channel(ChUsed(i)).Loc(:,1);
% end;
% figure
% pcntg = 0.5;
% for it = 1:2
%     %project 
%     CSpRSS = abs(UP'*Cfftp{it}*UP);
% %     a = (UP'*Cfft{it}*UP)./sqrt((diag(UP'*Cfft{it}*UP)*diag(UP'*Cfft{it}*UP)'));
%      a = (UP'*Cfft{it}*UP); %./sqrt((diag(UP'*Cfft{it}*UP)*diag(UP'*Cfft{it}*UP)'));
%     CSRSS = abs(imag(a));
%     
%     ind_max = find(CSpRSS(:)>pcntg*max(CSpRSS(:)));
%     subplot(2,3,(it-1)*3+1);
%     plot3(ChLoc(1,:),ChLoc(2,:),ChLoc(3,:),'.');
%     
%     hold on
%     for i=1:length(ind_max)
%         [ii jj]  = ind2sub(size(CSpRSS),ind_max(i));
%         plot3([ChLoc(1,ii) ChLoc(1,jj)],[ChLoc(2,ii) ChLoc(2,jj)],[ChLoc(3,ii) ChLoc(3,jj)],'Color','r');
%     end;
%     axis tight
%     axis off
%     ind_max = find(CSRSS(:)>0.8*max(CSRSS(:)));
%     subplot(2,3,(it-1)*3+2);
%     plot3(ChLoc(1,:),ChLoc(2,:),ChLoc(3,:),'.');
%     hold on
%     for i=1:length(ind_max)
%         [ii jj]  = ind2sub(size(CSRSS),ind_max(i));
%         plot3([ChLoc(1,ii) ChLoc(1,jj)],[ChLoc(2,ii) ChLoc(2,jj)],[ChLoc(3,ii) ChLoc(3,jj)],'Color','r');
%     end;
%     axis tight
%     axis off
%     
%     ind_max = find(PLI{it}(:)>0.8*max(PLI{it}(:)));
%     subplot(2,3,(it-1)*3+3);
%     plot3(ChLoc(1,:),ChLoc(2,:),ChLoc(3,:),'.');
%     
%     hold on
%     for i=1:length(ind_max)
%         [ii jj]  = ind2sub(size(PLV{it}),ind_max(i));
%         plot3([ChLoc(1,ii) ChLoc(1,jj)],[ChLoc(2,ii) ChLoc(2,jj)],[ChLoc(3,ii) ChLoc(3,jj)],'Color','r');
%     end;
%     axis tight
%     axis off
% 
% end;
return
CSp = real(CrossSpecTime{1}-Upwr{1}*(Upwr{1}'*CrossSpecTime{1}));
CS = real(CrossSpecTime{1});
CSN = real(CrossSpecTimeNoise{1});
CSD = real(CrossSpecTimeData{1});
CSpMax = max(abs(CSp),[],2);
CSMax = max(abs(CS),[],2);
indp = find(CSpMax>0.35*max(CSpMax));
ind = find(CSMax>0.35*max(CSMax));
indu = unique([ind;indp]);
figure
imagesc([CSN(indu,:) CS(indu,:) CSp(indu,:) abs(CSp(indu,:)-CSD(indu,:))]);

indp = find(max(abs(CSp),[],2))>0.5*max(abs(CSp(:)));
ind = find(max(abs(CS),[],2))>0.5*max(abs(CS(:)));


imagesc([real(CrossSpecTime{1}-Upwr{1}*(Upwr{1}'*CrossSpecTime{1})) real(CrossSpecTime{1})])


