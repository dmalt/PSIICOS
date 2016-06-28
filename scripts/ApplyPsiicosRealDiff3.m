% --------------------------- Set parameters --------------------------------------------- %
% ProtocolDir = 'C:/brainstorm_db/PSIICOS/data/';
ProtocolDir = '/home/dmalt/PSIICOS_osadtchii/data/';
% MEG_sensors_file = '/home/dmalt/ps/MEGSensors.mat';
MEG_sensors_file = './MEGSensors.mat';
% MEG_sensors_file = '/home/meg/osad/psiicos/MEGSensors.mat';
%  These flags are considered only when we load data from BST Protocol  %
%  and don't work if we load data from .mat file %
bUseHR = false;
bKeepLR = true;
bClearHM = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ChUsed = 1:306; ChUsed(3:3:end) = [];
TimeRange = [0, 0.700];
Conditions = {'1','2','4'}; % '2','4'};
Ncond = length(Conditions);
Band = [15 30];
BandName = 'beta';
%Band = [8 12];
bLoadTrials = true;
bComputePLI = false;
Fsamp = 500;
% [b,a] = butter(5, Band / (Fsamp / 2));
N_subjects = 10; % Need to write a function that will figure it out from the protocol
% ---------------------------------------------------------------------------------------- %
if exist(['./ConData_', BandName,'_CT.mat'], 'file')
    fprintf(['Loading ./ConData_', BandName, '_CT.mat; this might take a while ...\n']);
    load(['./ConData_', BandName,'_CT.mat']);
else
    if exist(['./ConData_', BandName, '.mat'], 'file')
        fprintf(['Loading ./ConData_', BandName, '.mat; this might take a while ...\n']);
        load(['./ConData_', BandName, '.mat']);
    else
        if exist('./10SubjData.mat', 'file')
            fprintf('Loading data from ./10SubjData.mat. This might take a while...\n')
            load('./10SubjData.mat');
        else
            % ---------------- Run brainstorm to read protocol info ----------------------------------- %
            brainstorm_path = '/home/dmalt/fif_matlab/brainstorm3/brainstorm';
            run( [brainstorm_path, '(''nogui'')'] );    % Start brainstorm without graphical interface
            % Protocol = bst_get('ProtocolStudies','PSIICOS');
            Protocol = bst_get('ProtocolStudies', 'PSIICOS_osadtchii');
            run( [brainstorm_path, '(''stop'')'] );    % Stop brainstorm 
            % ----------------------------------------------------------------------------------------- %

            clear ConData;
            fprintf('Loading real data from BST database.. \n');

            %---------- Load head models for subjects from brainstorm folders ------------------------------- %
            ConData = LoadHeadModels(Conditions, ProtocolDir, Protocol, bUseHR);

            % -------- Reduce tangent dimension and transform into virtual sensors -------------------------- %
            ConData = ReduceDimensions(ConData, ChUsed, bUseHR, bKeepLR, bClearHM);

            % --------------- Load trials from brainstorm ---------------------------------------%
            ConData = LoadTrials(ConData, Protocol, Conditions, bLoadTrials, ProtocolDir, ChUsed);

            disp('Saving ... \n');
            % save('c:\mywriteups\irAPMusicPaper\10SubjData.mat', '-v7.3');
            save('./10SubjData.mat', 'ConData', '-v7.3');
        end
        % ---------------- Band-pass filter and crop the data ------------------------ %
        ConData = BandPassFilter(ConData, Band, TimeRange, Fsamp);
        % ConData = ClearTrials(ConData);
        save(['./ConData_', BandName, '.mat'], 'ConData', '-v7.3');
    end
    % ---------- Compute cross-spectra for data and save them in .mat file ------------- %
    ConData = ComputeCrossSpectra(ConData);
    save(['./ConData_', BandName, '_CT.mat'], 'ConData', '-v7.3');
end

% -------- Load channel locations ------------------ %
ChLoc = ReadChannelLocations(MEG_sensors_file, ChUsed);
% ---- Split all cross-spectra by conditions ------- %
% In C we store cross-spectra we actually're going to
% use. C is a cell array of structures; for each 
% subject it stores cross-spectrum projected from volume
% conduction and from first condition, and matrix for
% going back to real sensors. 
[C1, C2, C4] = PrepCondCT(ConData, N_subjects);

TimeStep = (TimeRange(2) - TimeRange(1)) / (size(ConData{1}.CrossSpecTime, 2) - 1);
dtimes = TimeRange(1):TimeStep:TimeRange(2);
ind90 = TimeAsIndex(dtimes, 0.08);
ind110 = TimeAsIndex(dtimes, 0.14);

for s = 1:N_subjects
    % C{s}.CT = sum(C2{s}.CTfrom1, 2);
    % CT = sum(C2{s}.CT_Ind, 2);
    CT = sum(C1{s}.CT_Ind(:, 90:110),2);
    C{s}.CT = reshape(CT, sqrt(length(CT)), sqrt(length(CT)));
    C{s}.UP = C1{s}.UP;
    C{s}.G = C1{s}.G;
    C{s}.Upwr = C1{s}.Upwr;
end

RAPIts = 120;
Rnk = 350;
iSubj = 1;
[indep_topo, c_ss_hat, PVU, SubC, INDrap, Cp, Upwr] = RAP_PSIICOS_Fast(C{iSubj}.CT, C{iSubj}.G, RAPIts, Rnk, C{iSubj}.Upwr);

% Plot original cross-spectrum and residual after RAP_PSIICOS 
% side by side in similar colormap
CT = C1{iSubj}.CT_IndP;
[u s v] = svd(indep_topo, 0);
CT_res = CT - u * u' * CT;
CT = RestoreSensorDimension(CT, C{iSubj}.UP);
CT_res = RestoreSensorDimension(CT_res, C{iSubj}.UP);
drCT = abs(CT(:,:));
drCT_res = abs(CT_res(:,:));
bottom = min(min(drCT(:)), min(drCT_res(:) ));
top = max(max(drCT(:) ), max(drCT_res(:) ));
subplot(1,2,1);
imagesc(drCT,[bottom top]);
colorbar;
subplot(1,2,2);
imagesc(drCT_res, [bottom top])
colorbar;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pairs = PlotConnections(C, ChLoc, 'real');

CSa = zeros(10, 10);
for s1 = 1:10
    for s2 = 1:10
        if(s1 ~= s2)
            CSa(s1, s2) = ConnectivitySimilarity(Pairs{s1}, Pairs{s2}, ChLoc);
        end
    end;
end;

return
% ################################################################################################## %
% ################################################################################################## %
% ################################################################################################## %
% ################################################################################################## %
% ################################################################################################## %
sc = 1;
Trials1 = zeros(size(ConData{sc}.UP,2),size(ConData{sc}.Trials,2),size(ConData{sc}.Trials,2));
for tr = 1:size(ConData{sc}.Trials,3)
    Trials(:,:,tr) = ConData{sc}.UP'*ConDataBand{sc}.Trials(:,:,tr);
end;
Ep1 = squeeze(Trials(Pairs{1}(1,1),:,:))';
Ep2 = squeeze(Trials(Pairs{1}(1,2),:,:))';
Ep1Ind = Ep1-repmat(mean(Ep1,1),size(Ep1,1),1);
Ep2Ind = Ep2-repmat(mean(Ep1,1),size(Ep1,1),1);


for sc = 1:N_conditions_total
    SubjInd = fix((sc-1)/Ncond)*Ncond+1;
    ConData{sc}.CrossSpecTimeNoVC = ProjectAwayFromPower(ConData{sc}.CrossSpecTime,ConData{SubjInd}.G2dLRU);
end;


CT3  = ProjectAwayFromPower(ConData{3}.CrossSpecTime,G2dLRU);
[u2 s2 v2 ] = svd(CT2, 'econ');
CT3no2 = CT3 - u2(:,1:20) * (u2(:, 1:20)' * CT3);

[ Cs, CT, IND, Upwr] = RAPPSIICOSTime2Cond(ConData{3}.CrossSpecTime, ...
                 ConData{2}.CrossSpecTime, CT2, 20, G2dLRU ,1, 350, 5);





%ConData{2}.CrossSpecTime = CrossSpectralTimeseries( ConData{2}.Trials);
%ConData{3}.CrossSpecTime = CrossSpectralTimeseries( ConData{3}.Trials);
%[Qpsiicos, IND, CpProjs ] = RAPPSIICOS(ConData{2}.CrossSpec-ConData{1}.CrossSpec, G2dLRU,G2dHRU,4);
%[Qpsiicos, IND, CpProjs ] = RAPPSIICOS(ConData{3}.CrossSpec-ConData{1}.CrossSpec, G2dLRU,G2dHRU,3);
%[Q3vs1, IND, CpProjs3vs1, Upwr ] = PSIICOS(ConData{3}.CrossSpec-ConData{1}.CrossSpec, G2dLRU);
%[Q2vs1, IND, CpProjs2vs1, Upwr ] = PSIICOS(ConData{2}.CrossSpec-ConData{1}.CrossSpec, G2dLRU);
%[Q3, IND, CpProjs3, Upwr ] = PSIICOS(ConData{3}.CrossSpec, G2dLRU);
%[Q2, IND, CpProjs2, Upwr ] = PSIICOS(ConData{2}.CrossSpec, G2dLRU);
return

[Q3vs1T, IND, CpProjs3vs1, Upwr] = RAPPSIICOSTime(ConData{3}.CrossSpecTime - ConData{1}.CrossSpecTime, G2dLRU,4);
%[Cs3, Ps, INDdics] = iDICS(ConData{3}.CrossSpec, G2dLRU);
%[Cs1, Ps, INDdics] = iDICS(ConData{1}.CrossSpec, G2dLRU);

return
Cons = [1,3];
Ntr = size(ConData{3}.Trials,3);
for mc = 1:100
    trials = fix(0.99*rand(1,Ntr)*Ntr+1);
    CrossSpec = CrossSpectralMatrix(ConData{3}.Trials(:,:,trials),Band,500);
    [Qpsiicosmc, IND, CpProjs ] = RAPPSIICOS(CrossSpec-ConData{1}.CrossSpec, G2dLRU,G2dHRU,4);
    fname = sprintf('qpsiicos_mc_trial_%d.mat',mc);
    save(fname,'Qpsiicosmc');
%    QQmc{mc} = Qpsiicosmc;
end;
    
return
Ctx = load('D:\Brainstorm_db\PSIICOS\anat\0003_pran\tess_cortex_concat_2000V.mat');
%4175428 29
figure;
hctx  = trisurf(Ctx.Faces,Ctx.Vertices(:,1),Ctx.Vertices(:,2),Ctx.Vertices(:,3),'FaceColor',[0.1,0.51,1], 'EdgeColor','none','FaceAlpha', 0.3);
hold on;
camlight left; lighting phong
camlight right; 
hold on;

cols = ['r','g','m','y','k','c']
D =sum(Cs{1},2);
R = ConData{1}.HM_LR.GridLoc;
for k=1:size(D,2)
    Dmx = max(D(:,k));
    ind = find(D(:,k)>0.93*Dmx);
    for i=1:length(ind)
        h = line([R(IND(ind(i),1),1) R(IND(ind(i),2),1)],[R(IND(ind(i),1),2) R(IND(ind(i),2),2)],[R(IND(ind(i),1),3) R(IND(ind(i),2),3)] );
        plot3(R(IND(ind(i),1),1),R(IND(ind(i),1),2),R(IND(ind(i),1),3),[cols(k) '.']);
        plot3(R(IND(ind(i),2),1),R(IND(ind(i),2),2),R(IND(ind(i),2),3),[cols(k) '.']);
        set(h,'Color',cols(k),'LineWidth',2);
    end;
end;

clear CP;
CP = Cp{1}; 
[Qs key] = sort(QpsiicosP{1});
INDs = IND{1}(key,:);
tmp0 = CP;
a0 = norm(tmp0(:));
aa = zeros(1,100);
VV = [];
for r=1:100
    ii = INDs(end-k+1,1);
    jj = INDs(end-k+1,2);
    range_i = ii*2-1:ii*2;
    range_j = jj*2-1:jj*2;
    gi = G2dU(:,range_i);
    gj = G2dU(:,range_j);
    V = zeros(73^2, 4);
    Vre = zeros(73^2, 4);
    Vim = zeros(73^2, 4);
    k = 1;
    for i=1:2
        for j=1:2
            gg =bsxfun(@times,gi(:,i),gj(:,j)'); 
            v = gg+gg';
            Vre(:,k) = v(:);
            v = gg-gg';
            Vim(:,k) = v(:);
            k = k+1;
        end;
    end;
    VV= [VV V];
    [u s v] = svd(VV,'econ');
    c = u'*CP(:);
    CPp = reshape(CP(:)-u*c,73,73);
%    aare(k) = norm(real(tmp1(:)))/norm(real(tmp0(:)));
%    aaim(k) = norm(imag(tmp1(:)))/norm(imag(tmp0(:)));
    aa(r) = norm((CPp(:)))/norm((CP(:)));
end;

for mc=1:40
    fname = sprintf('qpsiicos_mc_trial_%d.mat',mc);
    h = load(fname);

    D = h.Qpsiicosmc;
    R = ConData{1}.HM_LR.GridLoc;
    for k=1:size(D,2)
        [Dmx ind] = max(D(:,k));
        for i=1:length(ind)
         h = line([R(IND(ind(i),1),1) R(IND(ind(i),2),1)],[R(IND(ind(i),1),2) R(IND(ind(i),2),2)],[R(IND(ind(i),1),3) R(IND(ind(i),2),3)] );
        plot3(R(IND(ind(i),1),1),R(IND(ind(i),1),2),R(IND(ind(i),1),3),[cols(k) '.']);
        plot3(R(IND(ind(i),2),1),R(IND(ind(i),2),2),R(IND(ind(i),2),3),[cols(k) '.']);
        set(h,'Color',cols(k));
        end;
    end;
    mc
end



