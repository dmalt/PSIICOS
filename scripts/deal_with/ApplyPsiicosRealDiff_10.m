%% Set params
% FolderName = 'C:/brainstorm_db/PSIICOS/data/';
FolderName = '/home/dmalt/PSIICOS_osadtchii/data/';
bUseHR = false;
ChUsed = 1:306; ChUsed(3:3:end) = [];
Conditions = {'1','2','4'}; % '2','4'};
Ncond = length(Conditions);
Band = [18 21];
%Band = [8 12];
bLoadTrials = true;
bComputePLI = false;
bStudyRank = false;
Fsamp = 500;
Protocol = bst_get('ProtocolStudies','PSIICOS');
clear ConData;
fprintf('Loading real data from BST database.. \n');
%% Load data and compute cross-spectral matrix 
ConditionsFound = 0;
clear ConData;
sc = 1;
for c = 1:length(Conditions)
    for s = 1:length(Protocol.Study)
        if(strcmp(Protocol.Study(s).Name,Conditions{c}))
            fprintf('Found study condition %s \n ', Conditions{c}); 
            for hm = 1:length(Protocol.Study(s).HeadModel)
                if(strcmp(Protocol.Study(s).HeadModel(hm).Comment,'Overlapping spheres_HR'))
                    ConData{sc}.HM_HR = load([FolderName Protocol.Study(s).HeadModel(hm).FileName]);
                else
                    ConData{sc}.HM_LR = load([FolderName Protocol.Study(s).HeadModel(hm).FileName]);
                end
            end;
            sc = sc+1;
        end;
    end;
end;
%% Reduce tangent dimension and transform into virtual sensors 
% the forward model is the same for both conditions
% so pick the first oneCOnData
GainSVDTh = 0.01;
Nch    = length(ChUsed);
for c = 1:length(ConData)
    ConData{c}.NsitesLR = size(ConData{c}.HM_LR.GridLoc,1);
    ConData{c}.G2dLR = zeros(Nch,ConData{c}.NsitesLR*2);
    % reduce tangent space
    range = 1:2;
    for i=1:ConData{c}.NsitesLR
        g = [ConData{c}.HM_LR.Gain(ChUsed,1+3*(i-1)) ...
             ConData{c}.HM_LR.Gain(ChUsed,2+3*(i-1)) ...
             ConData{c}.HM_LR.Gain(ChUsed,3+3*(i-1))];
        [u sv v] = svd(g);
        gt = g*v(:,1:2);
        ConData{c}.G2dLR(:,range) = gt*diag(1./sqrt(sum(gt.^2,1)));
        range = range + 2;
    end;
    
    %reduce sensor space
    [ug sg vg] = spm_svd(ConData{c}.G2dLR*ConData{c}.G2dLR',GainSVDTh);
    ConData{c}.UP = ug';
    ConData{c}.G2dLRU = ConData{c}.UP*ConData{c}.G2dLR;
    
    if(bUseHR)
        ConData{c}.NsitesHR = size(ConData{c}.HM_HR.GridLoc,1);
        ConData{c}.G2dHR = zeros(Nch,ConData{c}.NsitesHR*2);
        % reduce tangent space
        range = 1:2;
        for i=1:ConData{c}.NsitesHR
            g = [ConData{c}.HM_HR.Gain(ChUsed,1+3*(i-1)) ...
                 ConData{c}.HM_HR.Gain(ChUsed,2+3*(i-1)) ...
                 ConData{c}.HM_HR.Gain(ChUsed,3+3*(i-1))];
            [u sv v] = svd(g);
            gt = g*v(:,1:2);
            ConData{c}.G2dHR(:,range) = gt*diag(1./sqrt(sum(gt.^2,1)));
            range = range + 2;
        end;
    end;
    c;
end;

sc = 1;
ConditionsFound = 0;
for c = 1:length(Conditions)
    for s = 1:length(Protocol.Study)
        if(strcmp(Protocol.Study(s).Name,Conditions{c}))
            fprintf('Found study condition %s \n : ', Conditions{c}); 
            ConData{sc}.NumTrials = length(Protocol.Study(s).Data);
            if(bLoadTrials)
                fprintf('Loading Trials (Max %d) : ', ConData{sc}.NumTrials); 
    %            UP = ConData{fix((sc-1)/Ncond)*Ncond+1}.UP;
                UP = ConData{sc}.UP;
                for t = 1:ConData{sc}.NumTrials
                    aux = load([FolderName Protocol.Study(s).Data(t).FileName]);
                    if(t==1)
                         ConData{sc}.Trials = zeros(size(UP,1),length(aux.Time));
                         ConData{sc}.Time = aux.Time;
                         ConData{sc}.Fsamp = 1./(aux.Time(2)-aux.Time(1));
                    end;
                    ConData{sc}.Trials(:,:,t) = UP*aux.F(ChUsed,:);
                    if t>1
                        for tt=0:log10(t-1)
                            fprintf('\b'); % delete previous counter display
                        end
                    end
                    fprintf('%d', t);
                end; % trials t
                fprintf(' -> Done\n');
            end;
            sc = sc+1;
         end;
    end;
end;
% disp('Saving ... \n');
% save('c:\mywriteups\irAPMusicPaper\10SubjData.mat','-v7.3');
% return

% load('c:\mywriteups\irAPMusicPaper\10SubjData.mat');
%for sc = 1:10
%    ConData{sc} = [];
%end;


% do band-pass filtering and create ConDataBand
[bfir, afir] = fir1(128,Band/(Fsamp/2),'bandpass');

TimeRange(1) = -0.2;
TimeRange(2) = 0.6;
for sc = 1:length(ConData)
    for t = 1:size(ConData{sc}.Trials,3)
        [~, ind0] =min(abs(aux.Time-TimeRange(1)));
        [~, ind1] =min(abs(aux.Time-TimeRange(2)));
        T = ind1-ind0+1; 
        tmp = filtfilt(bfir,afir,(ConData{sc}.Trials(:,:,t))')';
        ConDataBand{sc}.Trials(:,:,t) = 1e12*tmp(:,ind0:ind1);
    end;
    sc
end;
bComputePLI = true;
for sc = 1:length(ConDataBand)
    fprintf('%d Computing cross-spectral matrix ....\n' , sc); 
    ConDataBand{sc}.CrossSpecTime = CrossSpectralTimeseries( ConDataBand{sc}.Trials);
    ConDataBand{sc}.CrossSpecTimeInd = CrossSpectralTimeseries( ConDataBand{sc}.Trials,true);
    % compute their projected versions
    [ConDataBand{sc}.CrossSpecTimeP, ConDataBand{sc}.Upwr] = ps.ProjectFromSlComplete(ConDataBand{sc}.CrossSpecTime, ConData{sc}.G2dLRU,350);
    ConDataBand{sc}.CrossSpecTimeIndP = ConDataBand{sc}.CrossSpecTimeInd - ConDataBand{sc}.Upwr*ConDataBand{sc}.Upwr'*ConDataBand{sc}.CrossSpecTimeInd;
    %UP
    if(bComputePLI)
        % UP is stored only for ConData and is the same as in ConDataBand,
        % so use the former
        Trials = zeros(size(ConData{sc}.UP,2),size(ConDataBand{sc}.Trials,2),size(ConDataBand{sc}.Trials,2));
        for tr = 1:size(ConData{sc}.Trials,3)
            Trials(:,:,tr) = ConData{sc}.UP'*ConDataBand{sc}.Trials(:,:,tr);
        end;
        ConDataBand{sc}.wPLI =  wPLIMatrix(Trials(:,1:256,:),Band,Fsamp,true);
    end;
end;
% study the amount of power reduction with VC projection rank increase
if(bStudyRank)
    RNK = [10 50 200 350 500];
    clear DS;
    %for rnk = 1:length(RNK)
    for rnk = 1:length(RNK)
        for sc = 1:length(ConDataBand)
            fprintf('%d Computing cross-spectral matrix ....\n' , sc); 
            % compute their projected versions
            [ConDataBand{sc}.CrossSpecTimeP, ConDataBand{sc}.Upwr,ds] = ps.ProjectFromSlComplete(ConDataBand{sc}.CrossSpecTime, ConData{sc}.G2dLRU,RNK(rnk));
            ConDataBand{sc}.CrossSpecTimeIndP = ConDataBand{sc}.CrossSpecTimeInd - ConDataBand{sc}.Upwr*ConDataBand{sc}.Upwr'*ConDataBand{sc}.CrossSpecTimeInd;
            RatInd(sc,rnk) = sum(abs(real(ConDataBand{sc}.CrossSpecTimeIndP(:))))/sum(abs(real(ConDataBand{sc}.CrossSpecTimeInd(:))));
            RatTot(sc,rnk) = sum(abs(real(ConDataBand{sc}.CrossSpecTimeP(:))))/sum(abs(real(ConDataBand{sc}.CrossSpecTime(:))));
            DS{rnk}(sc,:)=ds(:)';
        end
   end
end
UpperTri10 = [];
N = 10;
for i=1:N
    for j=i:N
        UpperTri10= [UpperTri10 sub2ind([N,N],i,j)];
    end;
end;

UpperTri8 = [];
N = 8;
for i=1:N
    for j=i:N
        UpperTri8= [UpperTri8 sub2ind([N,N],i,j)];
    end;
end;

load('c:\MyWriteups\iRAPMusicPaper\Simulations\MEGSensors.mat');
for i = 1:length(ChUsed)
    ChLoc(:,i) = MEGSensors.Channel(ChUsed(i)).Loc(:,1);
end;

for i=1:length(ChUsed)
     POS = MEGSensors.Channel(ChUsed(i)).Loc;
    ChOr(:,i) = (POS(:,2)-POS(:,3))/norm(POS(:,2)-POS(:,3));
end;


UpperTri = [];
N = 204;
for i=1:N
    for j=i:N
        UpperTri = [UpperTri sub2ind([N,N],i,j)];
    end;
end;

% calculate individual time-resolved cross-spectra averaged over all sensors 
S  = [1 2 3 6 7 8 9 10 4 5];
range = 1:50;
clear m4 m2 ss2 ss4 CS4vs1 CS2vs1 CS4x2 RNG;

for r = 1:20
    A2vs1 = zeros(size(ConDataBand{1}.CrossSpecTimeIndP));
    A4vs1 = zeros(size(ConDataBand{1}.CrossSpecTimeIndP));
    for s=1:10
        A2vs1 = zeros(size(ConDataBand{s}.CrossSpecTimeIndP));
        A4vs1 = zeros(size(ConDataBand{s}.CrossSpecTimeIndP));
        [u ss v] = svd(ConDataBand{s}.CrossSpecTimeIndP);
        A2vs1 =  ConDataBand{10+s}.CrossSpecTimeIndP - u(:,1:6)*u(:,1:6)'*ConDataBand{s}.CrossSpecTimeIndP;
        A4vs1 =  ConDataBand{20+s}.CrossSpecTimeIndP - u(:,1:6)*u(:,1:6)'*ConDataBand{s}.CrossSpecTimeIndP;
        rnk = size(ConData{s}.UP,1);
        CS_4vs1{s}  = ConData{s}.UP'*reshape(sum(A4vs1(:,range),2),rnk,rnk)*ConData{s}.UP;
        CS_2vs1{s}  = ConData{s}.UP'*reshape(sum(A2vs1(:,range),2),rnk,rnk)*ConData{s}.UP;
    end;
    for p=1:10
   
    pcntg = 0.002;
    for s = 1:10
        M = abs(real(CS_4vs1{s}));
        [aux, key_srt] = sort(M(:));
    %    ind_max = key_srt(fix((1-pcntg)*length(key_srt)):end);
        ind_max = key_srt(end-100*p:end);
        Pairs4vs1{s} = [];
        for i=1:length(ind_max)
            [ii jj]  = ind2sub(size(M),ind_max(i));
            Pairs4vs1{s} = [Pairs4vs1{s};[ii jj]];
        end;
        M = abs(real(CS_2vs1{s}));
        [aux, key_srt] = sort(M(:));
    %    ind_max = key_srt(fix((1-pcntg)*length(key_srt)):end);
        ind_max = key_srt(end-100*p:end);
        Pairs2vs1{s} = [];
        for i=1:length(ind_max)
            [ii jj]  = ind2sub(size(M),ind_max(i));
            Pairs2vs1{s} = [Pairs2vs1{s};[ii jj]];
        end;
    end

    CS4vs1{r} = zeros(10,10);
    CS2vs1{r} = zeros(10,10);
    CS4x2{r} = zeros(10,10);
    for s1 = 1:10
        for s2 = s1+1:10
            if(s1~=s2)
                CS4vs1{r}(s1,s2) = ConnectivitySimilarity(Pairs4vs1{S(s1)},Pairs4vs1{S(s2)},ChLoc);
                CS4x2{r}(s1,s2) = ConnectivitySimilarity(Pairs4vs1{S(s1)},Pairs2vs1{S(s2)},ChLoc);
                CS2vs1{r}(s1,s2) = ConnectivitySimilarity(Pairs2vs1{S(s1)},Pairs2vs1{S(s2)},ChLoc);
            end
        end;
    end;
    %RP2im(p) = mean(CS2vs1{r}(UpperTri10));
    %RP4im(p) = mean(CS4vs1{r}(UpperTri10));
    f = CS2vs1{r}(1:8,1:8);
    RP2re(p) = mean(f(UpperTri8));
    f = CS4vs1{r}(1:8,1:8);
    RP4re(p) = mean(f(UpperTri8));
    
    end
%     [pwi(r)]=ranksum(CS4vs1(UpperTri10),CS2vs1(UpperTri10));
%     [pac1(r)]=ranksum(CS4vs1(UpperTri10),CS4x2(UpperTri10));
%     [pac2(r)]=ranksum(CS2vs1(UpperTri10),CS4x2(UpperTri10));
    m4(r) = mean(CS4vs1{r}(UpperTri10));
    m2(r) = mean(CS2vs1{r}(UpperTri10));
    m42(r) = mean(CS4x2{r}(UpperTri10));
    
    RNG(r,:) =range;
    ss4(r) = std(CS4vs1{r}(UpperTri10));
    ss2(r) = std(CS2vs1{r}(UpperTri10));
    ss42(r) = std(CS4x2{r}(UpperTri10));
    
    range = range+20;
    r
end

[mxval,mxpos] = max(re.m2);
range = 201:250;

A2vs1 = zeros(size(ConDataBand{1}.CrossSpecTimeIndP));
A4vs1 = zeros(size(ConDataBand{1}.CrossSpecTimeIndP));
for s=1:10
    A2vs1 = zeros(size(ConDataBand{s}.CrossSpecTimeIndP));
    A4vs1 = zeros(size(ConDataBand{s}.CrossSpecTimeIndP));
    [u ss v] = svd(ConDataBand{s}.CrossSpecTimeIndP);
    A2vs1 =  ConDataBand{10+s}.CrossSpecTimeIndP - u(:,1:6)*u(:,1:6)'*ConDataBand{s}.CrossSpecTimeIndP;
    A4vs1 =  ConDataBand{20+s}.CrossSpecTimeIndP - u(:,1:6)*u(:,1:6)'*ConDataBand{s}.CrossSpecTimeIndP;
    rnk = size(ConData{s}.UP,1);
    CS_4vs1{s}  = ConData{s}.UP'*reshape(sum(A4vs1(:,range),2),rnk,rnk)*ConData{s}.UP;
    CS_2vs1{s}  = ConData{s}.UP'*reshape(sum(A2vs1(:,range),2),rnk,rnk)*ConData{s}.UP;
end;

for p=1:10

    pcntg = 0.002;
    for s = 1:10
        alpha = max(abs(imag(CS_4vs1{s}(:))))/max(abs(real(CS_4vs1{s}(:)))); 
        M = abs(real(CS_4vs1{s})*alpha +sqrt(-1)*imag(CS_4vs1{s}));
        [aux, key_srt] = sort(M(:));
    %    ind_max = key_srt(fix((1-pcntg)*length(key_srt)):end);
        ind_max = key_srt(end-100*p:end);
        Pairs4vs1{s} = [];
        for i=1:length(ind_max)
            [ii jj]  = ind2sub(size(M),ind_max(i));
            Pairs4vs1{s} = [Pairs4vs1{s};[ii jj]];
        end;
        %M = abs(real(CS_2vs1{s}));
        alpha = max(abs(imag(CS_2vs1{s}(:))))/max(abs(real(CS_2vs1{s}(:)))); 
        M = abs(real(CS_2vs1{s})*alpha +sqrt(-1)*imag(CS_2vs1{s}));

        [aux, key_srt] = sort(M(:));
    %    ind_max = key_srt(fix((1-pcntg)*length(key_srt)):end);
        ind_max = key_srt(end-100*p:end);
        Pairs2vs1{s} = [];
        for i=1:length(ind_max)
            [ii jj]  = ind2sub(size(M),ind_max(i));
            Pairs2vs1{s} = [Pairs2vs1{s};[ii jj]];
        end;
    end

        CS4vs1{r} = zeros(10,10);
        CS2vs1{r} = zeros(10,10);
        CS4x2{r} = zeros(10,10);
        for s1 = 1:10
            for s2 = s1+1:10
                if(s1~=s2)
                    CS4vs1{r}(s1,s2) = ConnectivitySimilarity(Pairs4vs1{S(s1)},Pairs4vs1{S(s2)},ChLoc);
                    CS4x2{r}(s1,s2) = ConnectivitySimilarity(Pairs4vs1{S(s1)},Pairs2vs1{S(s2)},ChLoc);
                    CS2vs1{r}(s1,s2) = ConnectivitySimilarity(Pairs2vs1{S(s1)},Pairs2vs1{S(s2)},ChLoc);
                end
            end;
        end;
        %RP2im(p) = mean(CS2vs1{r}(UpperTri10));
        %RP4im(p) = mean(CS4vs1{r}(UpperTri10));
        f = CS2vs1{r}(1:8,1:8);
        RP2abs(p) = mean(f(UpperTri8));
        f = CS4vs1{r}(1:8,1:8);
        RP4abs(p) = mean(f(UpperTri8));
end


 M = abs(real(CS_2vs1{s}));
[aux, key_srt] = sort(M(:));
        
ind_max = key_srt(end-100:end);
Pairs4vs1{s} = [];
CosAngRe = [];
k=1;
for i=1:length(ind_max)
    [ii jj]  = ind2sub(size(M),ind_max(i));
    if( norm(ChLoc(:,ii)-ChLoc(:,jj)) <0.35)
        CosAngRe(k) = ChOr(:,jj)'*ChOr(:,ii);
        k = k+1;
    end
end

        
        
s = 10;
M = abs(imag(CS_2vs1{s}));
[aux, key_srt] = sort(M(:));
%    ind_max = key_srt(fix((1-pcntg)*length(key_srt)):end);
ind_max = key_srt(end-399:end);
k = 1;
clear fi_im;
for i=1:length(ind_max)
    [ii jj]  = ind2sub(size(M),ind_max(i));
    if(norm(ChLoc(:,ii)-ChLoc(:,jj))>0.05)
        fi_im(k) = angle(CS_2vs1{s}(ii,jj));
        k = k+1;
    end;
end;

M = abs(real(CS_4vs1{s}));
[aux, key_srt] = sort(M(:));
%    ind_max = key_srt(fix((1-pcntg)*length(key_srt)):end);
ind_max = key_srt(end-399:end);
k=1;
clear fi_re;
for i=1:length(ind_max)
    [ii jj]  = ind2sub(size(M),ind_max(i));
     if(norm(ChLoc(:,ii)-ChLoc(:,jj))>0.05)
        fi_re(k) = angle(CS_2vs1{s}(ii,jj));
        k = k+1;
    end;
end;

[him,aim] = hist(fi_im);
[hre,are] = hist(fi_re);
figure
plot(aim,him,'r');
hold on
plot(are,hre,'b');
legend('phase distribution of max(abs(imag))','phase distribution of max(abs(real))');




    
figure
imagesc(CSa)
UpperTri10 = [];
N = 10;
for i=1:N
    for j=i+1:N
        UpperTri10 = [UpperTri10 sub2ind([N,N],i,j)];
    end;
end;
for i=1:204
    for j=1:204
        dst(i,j) = norm(ChLoc(:,i)-ChLoc(:,j));
    end;
end;



    
%     A1P = (ConDataBand{10+s}.CrossSpecTimeIndP); % - ConDataBand{s}.CrossSpecTimeInd);
%     A2P = (ConDataBand{20+s}.CrossSpecTimeIndP);% - ConDataBand{s}.CrossSpecTimeInd);
%     A1 = (ConDataBand{10+s}.CrossSpecTimeInd); % - ConDataBand{s}.CrossSpecTimeInd);
%     A2 = (ConDataBand{20+s}.CrossSpecTimeInd);% - ConDataBand{s}.CrossSpecTimeInd);
%     for t = 1:T
%         N = size(ConData{10+s}.UP,1);
%         Acc21P{s}(:,:,t) = ConData{10+s}.UP'*reshape(A2P(:,t)-A1P(:,t),N,N)*ConData{10+s}.UP;
%         Acc21{s}(:,:,t) = ConData{10+s}.UP'*reshape(A2(:,t)-A1(:,t),N,N)*ConData{10+s}.UP;
%         Acc21VC{s}(:,:,t) =  ConData{20+s}.UP'*reshape((A2(:,t) - A2P(:,t)) -  (A1(:,t) - A1P(:,t)),N,N)*ConData{20+s}.UP;
%     end;
% end;

% close all
% figure
% for s = 1:9
%     subplot(3,3,s)
%     imagesc(sum(abs(real(Acc21{s})),3))
%     axis off
% end;
% figure
% for s = 1:9
%     subplot(3,3,s)
%     imagesc(sum(abs(real(Acc21P{s})),3))
%     axis off
% end;
% ACC = zeros(204,204);
% ACCP = zeros(204,204);
% COR = zeros(204,204);
% for s = 1:9
%     tmp = sum(abs(real(Acc21P{s})),3);
%     tmp  = tmp/sum(tmp(:));
%     ACCP = ACCP+tmp;
%     tmp = sum(abs(real(Acc21{s})),3);
%     tmp  = tmp/sum(tmp(:));
%     ACC = ACC+tmp;
%     COR = COR + sum(real(Acc21VC{s}(:,:,:)).*real(Acc21P{s}(:,:,:)),3)./(sqrt(sum(real(Acc21VC{s}(:,:,:)).*real(Acc21VC{s}(:,:,:)),3)*sum(real(Acc21P{s}(:,:,:)).*real(Acc21P{s}(:,:,:)),3)));
% end;
   
    
% ACCI = zeros(204,204);

% for s = 1:9
%     tmp = sum(abs(imag(Acc21{s})),3);
%     tmp  = tmp/sum(tmp(:));
%     ACCI = ACCI+tmp;
% end
% figure
% %Analysis of power independent entries of the cross-spectral matrix
% %Yellow — only significant real part present (detectable power independent zero phase locking)
% %Blue — only significant imaginary part present ( can't say about the phase on the near-diagonal locations as the real part is eliminated by the projection)
% %Dark red — both parts are present — non-pi/2 and non zero phase locking
% imagesc((abs(COR)<0.2).*(ACCP>4.5e-4)-(abs(COR)<0.2).*(ACCI>4.5e-4)+2*((abs(COR)<0.2).*(ACCP>4.5e-4)==(abs(COR)<0.2).*(ACCI>4.5e-4)).*(((abs(COR)<0.2).*(ACCP>4.5e-4)>0 ).*((abs(COR)<0.2).*(ACCI>4.5e-4))>0))


% for s=1:10
%     Acc21(s,:) = Acc21(s,:)/max(abs(Acc21(s,:)));
%     Acc2(s,:) = Acc2(s,:)/max(abs(Acc2(s,:)));
%     Acc1(s,:) = Acc1(s,:)/max(abs(Acc1(s,:)));
% end;

    
% figure
% plot(mean(abs(Acc2)),'r');
% hold on
% plot(mean(abs(Acc1)),'g');
% plot(mean(abs(Acc21)),'b');
% clear UP;
% for sc  = 1:10
%     [ug sg vg] = spm_svd(ConData{sc}.G2dLR*ConData{sc}.G2dLR',0.001);
%     UP{sc} = ug';
%     ProjRank(sc) = size(UP{sc},1);
% end;
% minSz = min(ProjRank);
% for sc  = 1:10
%     Ulf{sc} = UP{sc}(1:minSz,:);
%     G2dLRU{sc}= Ulf{sc}*ConData{sc}.G2dLR;
% end;

% %study the effect if systematic error
% for sc = 1:10
%     [Upwr{sc},ES{sc},A{sc}] = ProjectorOnlyAwayFromPowerComplete(G2dLRU{sc},350);
% end

% for s1 = 1:10
%     for s2 = 1:10
%         CTVCp  = A{s1}-Upwr{s2}*(Upwr{s2}'*A{s1});        
%         rat(s1,s2) = norm(CTVCp(:))/norm(A{s1}(:));
%     end;
%     s1
% end;
% figure
% imagesc(rat); colorbar
% % study random error
% sc = 1;
% NoiseLevel = [0 0.01 0.1 0.15 0.25 0.5 0.75 1 ];
% for sc = 1:10
%     G2dLRU{sc} = G2dLRU{sc}(:,1:4000);
% end;
% NoiseStr = zeros(size(G2dLRU{sc}));
% for sc = 2:3
%     NoiseStr = NoiseStr+G2dLRU{sc};
% end;
% NoiseStr = G2dLRU{3}-G2dLRU{2};
% NoiseStr = (NoiseStr - mean(NoiseStr(:)))/std(NoiseStr(:));

% [hst,argst] = hist(NoiseStr(:),30);
% x = argst;
% y = cumsum(hst)/sum(hst);
% cs = spline(y,linspace(0,1,30));
% NoiseRND = ppval(cs,rand(size(G2dLRU{sc})));
% NoiseRND = (NoiseRND - mean(NoiseRND(:)))/std(NoiseRND(:));
% [hst,argst] = hist(NoiseStr(:),30);
% [hrnd,argrnd] = hist(NoiseRND(:),30);
% figure
% plot(argst,hst);
% hold on
% plot(argrnd,hrnd,'r');


% clear ratRND  ratSTR;
% sc = 1;
% [~,~,A0{sc}] = ProjectorOnlyAwayFromPowerComplete(G2dLRU{sc},350);
% for nl = 1:length(NoiseLevel)
%     alpha = NoiseLevel(nl)*mean(abs(G2dLRU{sc}(:)));
%     [UpwrRND{sc},ESRND{sc},ARND{sc}] = ProjectorOnlyAwayFromPowerComplete(G2dLRU{sc}+alpha*NoiseRND,350);
%     CTVCp  = A0{sc}-UpwrRND{sc}*(UpwrRND{sc}'*A0{sc});        
%     ratRND(nl) = norm(CTVCp(:))/norm(A0{sc}(:));
%     [AttRe, AttIm] = MeasureAttenuationReIm( UpwrRND{sc},G2dLRU{sc},300);
%     AttReRND(nl) = AttRe;
%     AttImRND(nl) = AttIm;
   
%     [UpwrSTR{sc},ESSTR{sc},ASTR{sc}] = ProjectorOnlyAwayFromPowerComplete(G2dLRU{sc}+alpha*NoiseStr,350);
%     CTVCp  = A0{sc}-UpwrSTR{sc}*(UpwrSTR{sc}'*A0{sc});        
%     ratSTR(nl) = norm(CTVCp(:))/norm(A0{sc}(:));
%     [AttRe, AttIm] = MeasureAttenuationReIm( UpwrSTR{sc},G2dLRU{sc},300);
%     AttReSTR(nl) = AttRe;
%     AttImSTR(nl) = AttIm;
% end;

% figure
% subplot(2,1,1)
% plot(NoiseLevel,ratSTR,'ko-','LineWidth',1.5);
% hold on
% plot(NoiseLevel,ratRND,'k>-','LineWidth',1.5);
% grid
% xlabel('Noise level' );
% ylabel('ratio');
% title('VC attenuation');
% legend('Structured','Random')
% subplot(2,1,2)
% plot(NoiseLevel,AttReSTR,'ko-','LineWidth',1.5);
% hold on
% plot(NoiseLevel,AttReRND,'k>-','LineWidth',1.5);
% legend('Structured','Random')
% grid
% xlabel('Noise level' );
% ylabel('ratio');
% title('Re part of interaction component attenuation');





% CTImp  = imag(CT)-u(:,1:rnk)*(u(:,1:rnk)'*imag(CT));
%     CTRep  = real(CT)-u(:,1:rnk)*(u(:,1:rnk)'*real(CT));
%     CTVCp  = A-u(:,1:rnk)*(u(:,1:rnk)'*A);
%     nrmim(k) = norm(CTImp(:))/norm(imag(CT(:)));
%     nrmre(k) = norm(CTRep(:))/norm(real(CT(:)));
%     nrmvc(k) = norm(CTVCp(:))/norm(A(:));


% load('c:\MyWriteups\iRAPMusicPaper\Simulations\MEGSensors.mat');
% for i = 1:length(ChUsed)
%     ChLoc(:,i) = MEGSensors.Channel(ChUsed(i)).Loc(:,1);
% end;

% range = 75:150;
% figure
% pcntg = 2*1e-3;
% for s=1:10
%      C1 = ConDataBand{10+s}.CrossSpecTimeIndP(:,range)-ConDataBand{s}.CrossSpecTimeIndP(:,range);
%      C2 = ConDataBand{20+s}.CrossSpecTimeIndP(:,range)-ConDataBand{s}.CrossSpecTimeIndP(:,range);
%      [u ss2 v] = svd([real(C2) imag(C2)]);
%      C1but2 = C1-u(:,1:6)*u(:,1:6)'*C1;
%      [u ss1 v] = svd([real(C1) imag(C1)]);
%      C2but1 = C2-u(:,1:6)*u(:,1:6)'*C2;
%      C = sum(C1but2(:,1:50),2);

% %     C2 = ConDataBand{10+s}.CrossSpecTimeIndP(:,75:200);
% %     C1 = ConDataBand{s}.CrossSpecTimeIndP(:,75:200);
% %     
% %     [u ss v] = svd([real(C1) imag(C1)]);
% %     C2but1 = C2-u(:,1:15)*u(:,1:15)'*C2;
% %     C = sum(C2but1(:,1:50),2);
% %    
%     Csq = reshape(C,size(ConData{10+s}.UP,1),size(ConData{10+s}.UP,1));
    
%     D21{s} = ConData{10+s}.UP'*Csq*ConData{10+s}.UP;
%     M = abs(real(D21{s}));
%     %M = (ConDataBand{20+s}.wPLI-ConDataBand{s}.wPLI)-(ConDataBand{10+s}.wPLI-ConDataBand{s}.wPLI);
%     [aux, key_srt] = sort(M(:));
%     ind_max = key_srt(fix((1-pcntg)*length(key_srt)):end);
%     th = aux(fix((1-pcntg)*length(key_srt)));

%     h = subplot(2,5,s)
%       plot3(ChLoc(1,:),ChLoc(2,:),ChLoc(3,:),'.');

%       hold on
%       Pairs{s} = [];
%        for i=1:length(ind_max)
%           [ii jj]  = ind2sub(size(D21{s}),ind_max(i));
%           Pairs{s} = [Pairs{s};[ii jj]];
%           plot3([ChLoc(1,ii) ChLoc(1,jj)],[ChLoc(2,ii) ChLoc(2,jj)],[ChLoc(3,ii) ChLoc(3,jj)],'Color','r');
%         end;
%         set(h,'View',[0 90])
%         axis tight
%         axis off
% end
% CSa = zeros(10,10);
% for s1 = 1:10
%     for s2 = 1:10
%         if(s1~=s2)
%             CSa(s1,s2) = ConnectivitySimilarity(Pairs{s1},Pairs{s2},ChLoc);
%         end
%     end;
% end;

% return
% sc = 1;
% Trials1 = zeros(size(ConData{sc}.UP,2),size(ConData{sc}.Trials,2),size(ConData{sc}.Trials,2));
% for tr = 1:size(ConData{sc}.Trials,3)
%     Trials(:,:,tr) = ConData{sc}.UP'*ConDataBand{sc}.Trials(:,:,tr);
% end;
% Ep1 = squeeze(Trials(Pairs{1}(1,1),:,:))';
% Ep2 = squeeze(Trials(Pairs{1}(1,2),:,:))';
% Ep1Ind = Ep1-repmat(mean(Ep1,1),size(Ep1,1),1);
% Ep2Ind = Ep2-repmat(mean(Ep1,1),size(Ep1,1),1);


% for sc = 1:length(ConData)
%     SubjInd = fix((sc-1)/Ncond)*Ncond+1;
%     ConData{sc}.CrossSpecTimeNoVC = ProjectAwayFromPower(ConData{sc}.CrossSpecTime,ConData{SubjInd}.G2dLRU);
% end;


% CT3  = ProjectAwayFromPower(ConData{3}.CrossSpecTime,G2dLRU);
% [u2 s2 v2 ] = svd(CT2,'econ');
% CT3no2 = CT3-u2(:,1:20)*(u2(:,1:20)'*CT3);

% [ Cs, CT, IND, Upwr] = RAPPSIICOSTime2Cond(ConData{3}.CrossSpecTime, ConData{2}.CrossSpecTime,CT2,20,G2dLRU ,1,350, 5);





% %ConData{2}.CrossSpecTime = CrossSpectralTimeseries( ConData{2}.Trials);
% %ConData{3}.CrossSpecTime = CrossSpectralTimeseries( ConData{3}.Trials);
% %[Qpsiicos, IND, CpProjs ] = RAPPSIICOS(ConData{2}.CrossSpec-ConData{1}.CrossSpec, G2dLRU,G2dHRU,4);
% %[Qpsiicos, IND, CpProjs ] = RAPPSIICOS(ConData{3}.CrossSpec-ConData{1}.CrossSpec, G2dLRU,G2dHRU,3);
% %[Q3vs1, IND, CpProjs3vs1, Upwr ] = PSIICOS(ConData{3}.CrossSpec-ConData{1}.CrossSpec, G2dLRU);
% %[Q2vs1, IND, CpProjs2vs1, Upwr ] = PSIICOS(ConData{2}.CrossSpec-ConData{1}.CrossSpec, G2dLRU);
% %[Q3, IND, CpProjs3, Upwr ] = PSIICOS(ConData{3}.CrossSpec, G2dLRU);
% %[Q2, IND, CpProjs2, Upwr ] = PSIICOS(ConData{2}.CrossSpec, G2dLRU);
% return

% [Q3vs1T, IND, CpProjs3vs1, Upwr ] = RAPPSIICOSTime(ConData{3}.CrossSpecTime-ConData{1}.CrossSpecTime, G2dLRU,4);
% %[Cs3, Ps, INDdics] = iDICS(ConData{3}.CrossSpec, G2dLRU);
% %[Cs1, Ps, INDdics] = iDICS(ConData{1}.CrossSpec, G2dLRU);

% return
% Cons = [1,3];
% Ntr = size(ConData{3}.Trials,3);
% for mc = 1:100
%     trials = fix(0.99*rand(1,Ntr)*Ntr+1);
%     CrossSpec = CrossSpectralMatrix(ConData{3}.Trials(:,:,trials),Band,500);
%     [Qpsiicosmc, IND, CpProjs ] = RAPPSIICOS(CrossSpec-ConData{1}.CrossSpec, G2dLRU,G2dHRU,4);
%     fname = sprintf('qpsiicos_mc_trial_%d.mat',mc);
%     save(fname,'Qpsiicosmc');
% %    QQmc{mc} = Qpsiicosmc;
% end;
    
% return
% Ctx = load('D:\Brainstorm_db\PSIICOS\anat\0003_pran\tess_cortex_concat_2000V.mat');
% %4175428 29
% figure;
% hctx  = trisurf(Ctx.Faces,Ctx.Vertices(:,1),Ctx.Vertices(:,2),Ctx.Vertices(:,3),'FaceColor',[0.1,0.51,1], 'EdgeColor','none','FaceAlpha', 0.3);
% hold on;
% camlight left; lighting phong
% camlight right; 
% hold on;

% cols = ['r','g','m','y','k','c']
% D =sum(Cs{1},2);
% R = ConData{1}.HM_LR.GridLoc;
% for k=1:size(D,2)
%     Dmx = max(D(:,k));
%     ind = find(D(:,k)>0.93*Dmx);
%     for i=1:length(ind)
%         h = line([R(IND(ind(i),1),1) R(IND(ind(i),2),1)],[R(IND(ind(i),1),2) R(IND(ind(i),2),2)],[R(IND(ind(i),1),3) R(IND(ind(i),2),3)] );
%         plot3(R(IND(ind(i),1),1),R(IND(ind(i),1),2),R(IND(ind(i),1),3),[cols(k) '.']);
%         plot3(R(IND(ind(i),2),1),R(IND(ind(i),2),2),R(IND(ind(i),2),3),[cols(k) '.']);
%         set(h,'Color',cols(k),'LineWidth',2);
%     end;
% end;

% clear CP;
% CP = Cp{1}; 
% [Qs key] = sort(QpsiicosP{1});
% INDs = IND{1}(key,:);
% tmp0 = CP;
% a0 = norm(tmp0(:));
% aa = zeros(1,100);
% VV = [];
% for r=1:100
%     ii = INDs(end-k+1,1);
%     jj = INDs(end-k+1,2);
%     range_i = ii*2-1:ii*2;
%     range_j = jj*2-1:jj*2;
%     gi = G2dU(:,range_i);
%     gj = G2dU(:,range_j);
%     V = zeros(73^2, 4);
%     Vre = zeros(73^2, 4);
%     Vim = zeros(73^2, 4);
%     k = 1;
%     for i=1:2
%         for j=1:2
%             gg =bsxfun(@times,gi(:,i),gj(:,j)'); 
%             v = gg+gg';
%             Vre(:,k) = v(:);
%             v = gg-gg';
%             Vim(:,k) = v(:);
%             k = k+1;
%         end;
%     end;
%     VV= [VV V];
%     [u s v] = svd(VV,'econ');
%     c = u'*CP(:);
%     CPp = reshape(CP(:)-u*c,73,73);
% %    aare(k) = norm(real(tmp1(:)))/norm(real(tmp0(:)));
% %    aaim(k) = norm(imag(tmp1(:)))/norm(imag(tmp0(:)));
%     aa(r) = norm((CPp(:)))/norm((CP(:)));
% end;

% for mc=1:40
%     fname = sprintf('qpsiicos_mc_trial_%d.mat',mc);
%     h = load(fname);

%     D = h.Qpsiicosmc;
%     R = ConData{1}.HM_LR.GridLoc;
%     for k=1:size(D,2)
%         [Dmx ind] = max(D(:,k));
%         for i=1:length(ind)
%          h = line([R(IND(ind(i),1),1) R(IND(ind(i),2),1)],[R(IND(ind(i),1),2) R(IND(ind(i),2),2)],[R(IND(ind(i),1),3) R(IND(ind(i),2),3)] );
%         plot3(R(IND(ind(i),1),1),R(IND(ind(i),1),2),R(IND(ind(i),1),3),[cols(k) '.']);
%         plot3(R(IND(ind(i),2),1),R(IND(ind(i),2),2),R(IND(ind(i),2),3),[cols(k) '.']);
%         set(h,'Color',cols(k));
%         end;
%     end;
%     mc
% end



