close all;
clear all;
LoadCovariance = 0;
Nmc = 100;
PermutePreVsPost = 1;
Left1Right2 = 2;

Fa0 = 8; %Hz;
Fa1 = 12; %Hz;
Fb0 = 16;
Fb1 = 24; % Hz;
Fg0 = 50; %Hz;
Fg1 = 80; % Hz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set pair of bands for this analysis 
F_bnd1_0 = Fb0;
F_bnd1_1 = Fb1;

F_bnd2_0 = Fg0;
F_bnd2_1 = Fg1;


cd C:\LocalData\Biomag2010Competition\dataset2preprocessed;

G0 = ReadEMSEMatrix('G_IS.bin','double');
CM = ReadEMSEMatrix('CM_IS.bin','float');
G = CM*G0;
GridPoints = ReadEMSEMatrix('Gridpoints_IS.bin','float');

NDatasets = 5;
DataFileName = cell(5);
if(Left1Right2==2)
    DataFileName{1} = 'dataset2_2_1_preproc_v6.mat';
    DataFileName{2} = 'dataset2_2_2_preproc_v6.mat';
    DataFileName{3} = 'dataset2_2_3_preproc_v6.mat';
    DataFileName{4} = 'dataset2_2_4_preproc_v6.mat';
    DataFileName{5} = 'dataset2_2_5_preproc_v6.mat';
else
    DataFileName{1} = 'dataset2_1_1_preproc_v6.mat';
    DataFileName{2} = 'dataset2_1_2_preproc_v6.mat';
    DataFileName{3} = 'dataset2_1_3_preproc_v6.mat';
    DataFileName{4} = 'dataset2_1_4_preproc_v6.mat';
    DataFileName{5} = 'dataset2_1_5_preproc_v6.mat';
end;

load(DataFileName{1});


Fs = data.fsample;
%set bands

% start analysis
Nepochs = length(data.trial);
Nch = size(data.trial{1},1);

C1_pre = zeros(Nch,Nch);
C2_pre = zeros(Nch,Nch);

C1_post = zeros(Nch,Nch);
C2_post = zeros(Nch,Nch);

Cbeta_pre = zeros(Nch,Nch);
Cbeta_post = zeros(Nch,Nch);

Npre    = 0;
Npost   = 0;
Nbeta_pre = 0;
Nbeta_post = 0;

% if requested will read covariances from this file
filename_cov = sprintf('CovMatr_%dto%dvs%dto%d_R_FFT.mat',F_bnd1_0,F_bnd1_1,F_bnd2_0,F_bnd2_1); 

R1UpdatePreBnd1 = [];
R1UpdatePreBnd2 = [];
if(LoadCovariance == 0)
    
    for ds = 1:NDatasets
        disp('New dataset');
        if(ds>1)
            clear data;
            load(DataFileName{ds});
            Nepochs = length(data.trial);
        end;
        for tr =1:Nepochs
            range_pre   = find((data.time{tr}<=-0.2)&(data.time{tr}>=-0.6));
            range_post  = find((data.time{tr}>1.1)&(data.time{tr}<1.5));

            Xfft = fft(data.trial{tr}(:,range_pre),[],2)/length(range_pre);
            F = data.fsample/length(range_pre-1)*(0:(length(range_pre)-1));
            range_bnd1 = find((F>=F_bnd1_0) & (F<=F_bnd1_1));
            range_bnd2 = find((F>=F_bnd2_0) & (F<=F_bnd2_1));
            if(~isempty(range_bnd1) && ~isempty(range_bnd2))
                XfftMeanBnd1 = mean(Xfft(:,range_bnd1),2);
                XfftMeanBnd2 = mean(Xfft(:,range_bnd2),2);
            
                C1_pre = C1_pre+XfftMeanBnd1*XfftMeanBnd1'+XfftMeanBnd2*XfftMeanBnd2';
                C2_pre = C2_pre+(XfftMeanBnd1+XfftMeanBnd2)*(XfftMeanBnd1+XfftMeanBnd2)';
                Npre    = Npre+1;
                
                R1UpdatePreBnd1(:,Npre) = XfftMeanBnd1; 
                R1UpdatePreBnd2(:,Npre) = XfftMeanBnd2; 
            end;
            
            
            Xfft = fft(data.trial{tr}(:,range_post),[],2)/length(range_post);
            F = data.fsample/length(range_post-1)*(0:(length(range_post)-1));
            range_bnd1 = find((F>=F_bnd1_0) & (F<=F_bnd1_1));
            range_bnd2 = find((F>=F_bnd2_0) & (F<=F_bnd2_1));
            
            if(~isempty(range_bnd1) && ~isempty(range_bnd2))
                XfftMeanBnd1 = mean(Xfft(:,range_bnd1),2);
                XfftMeanBnd2 = mean(Xfft(:,range_bnd2),2);
            
                C1_post = C1_post+XfftMeanBnd1*XfftMeanBnd1'+XfftMeanBnd2*XfftMeanBnd2';
                C2_post = C2_post+(XfftMeanBnd1+XfftMeanBnd2)*(XfftMeanBnd1+XfftMeanBnd2)';
                Npost    = Npost+1;
            
                R1UpdatePostBnd1(:,Npost) = XfftMeanBnd1; 
                R1UpdatePostBnd2(:,Npost) = XfftMeanBnd2; 
            end;
            
            % comupte betta covariance
            range_pre   = find((data.time{tr}<=-0.1)&(data.time{tr}>=-0.6));
            Xfft = fft(data.trial{tr}(:,range_pre),[],2)/length(range_pre);
            F = data.fsample/length(range_pre-1)*(0:(length(range_pre)-1));
            range_beta = find((F>=Fb0 )& (F<=Fb1));
            if(~isempty(range_beta))
               XfftMeanBeta = mean(Xfft(:,range_beta),2);
                Cbeta_pre = Cbeta_pre+XfftMeanBeta*XfftMeanBeta'+XfftMeanBeta*XfftMeanBeta';
                Nbeta_pre    = Nbeta_pre+1;
            end;
            
            range_post  = find((data.time{tr}>=0.4)& (data.time{tr}<=0.9));
            Xfft = fft(data.trial{tr}(:,range_post),[],2)/length(range_post);
            F = data.fsample/length(range_post-1)*(0:(length(range_post)-1));
            range_beta = find((F>=Fb0 )& (F<=Fb1));
            if(~isempty(range_beta))
               XfftMeanBeta = mean(Xfft(:,range_beta),2);
                Cbeta_post = Cbeta_post+XfftMeanBeta*XfftMeanBeta'+XfftMeanBeta*XfftMeanBeta';
                Nbeta_post    = Nbeta_post+1;
            end;
            tr
        end;
    end;
    
    C1_pre = C1_pre/Npre;
    C2_pre = C2_pre/Npre;
    C1_post = C1_post/Npost;
    C2_post = C2_post/Npost;
    Cbeta_pre = Cbeta_pre/Nbeta_pre;
    Cbeta_post = Cbeta_post/Nbeta_post;
    
    save(filename_cov, 'C1_pre', 'C2_pre', 'C1_post', 'C2_post', 'Cbeta_pre','Cbeta_post');
else
    load(filename_cov);
end;
%return;    

% covariance data are ready, 
% build beamformers
alpha = 0.01; % 0.01 is default for tihinv
iC1_pre =  tihinv(C1_pre);
iC2_pre =  tihinv(C2_pre);
iC1_post = tihinv(C1_post);
iC2_post = tihinv(C2_post);
iCbeta_post = tihinv(Cbeta_post);
iCbeta_pre = tihinv(Cbeta_pre);

Nsites   = size(G,2)/3;
Nch = size(G,1);

W1_pre   = cell(Nsites);
W2_pre   = cell(Nsites);
W1_post  = cell(Nsites);
W2_post  = cell(Nsites);
Wbeta_post  = cell(Nsites);
Wbeta_pre   = cell(Nsites);

for i = 1:Nsites
    W1_pre{i} = zeros(3,Nch);
    W2_pre{i} = zeros(3,Nch);
    W1_post{i} = zeros(3,Nch);
    W2_post{i} = zeros(3,Nch);
    Wbeta_post{i} = zeros(3,Nch);
    Wbeta_pre{i} = zeros(3,Nch);
end

range = 1:3;
for i = 1:Nsites
    denum = tihinv(G(:,range)'*iC1_pre*G(:,range));
    num = G(:,range)'*iC1_pre;
    W1_pre{i} = denum*num;
  
    denum = tihinv(G(:,range)'*iC1_post*G(:,range));
    num = G(:,range)'*iC1_post;
    W1_post{i} = denum*num;

    denum = tihinv(G(:,range)'*iC2_pre*G(:,range));
    num = G(:,range)'*iC2_pre;
    W2_pre{i} = denum*num;
  
    denum = tihinv(G(:,range)'*iC2_post*G(:,range));
    num = G(:,range)'*iC2_post;
    W2_post{i} = denum*num;
    
    denum = tihinv(G(:,range)'*iCbeta_pre*G(:,range));
    num = G(:,range)'*iCbeta_pre;
    Wbeta_pre{i} = denum*num;
    denum = tihinv(G(:,range)'*iCbeta_post*G(:,range));
    num = G(:,range)'*iCbeta_post;
    Wbeta_post{i} = denum*num;
      
    range = range+3;
    i
end;
% compute auto terms
Q = zeros(1,Nsites);
Qbetadesync = zeros(1,Nsites);
for i = 1:Nsites
    A = W1_pre{i}*C1_pre* W1_pre{i}';
    [~, s ~] = svd(A);
    S1_pre = s(1,1);

    A = W1_post{i}*C1_post* W1_post{i}';
    [~, s, ~] = svd(A);
    S1_post = s(1,1);

    A = W2_pre{i}*C2_pre* W2_pre{i}';
    [~ ,s, ~] = svd(A);
    S2_pre = s(1,1);

    A = W2_post{i}*C2_post* W2_post{i}';
    [~, s, ~] = svd(A);
    S2_post = s(1,1);
    Q(i) = (S2_post/S1_post)/(S2_pre/S1_pre);
    
    A = Wbeta_pre{i}*Cbeta_pre* Wbeta_pre{i}';
    [~, s ~] = svd(A);
    Sbeta_pre = s(1,1);

    A = Wbeta_post{i}*Cbeta_post* Wbeta_post{i}';
    [~, s ~] = svd(A);
    Sbeta_post = s(1,1);
    Qbetadesync(i) = Sbeta_pre/Sbeta_post;
end;

Sites = GridPoints(1:3:end,1:3);

Z = unique(Sites(:,3));
Nlayers = length(Z);
M = 0.5*ones(17,17,Nlayers);

close all;
figure
colorbar
K = 150;
M = 0.5*ones(20,20*Nlayers);    
IND = [];
for i=1:Nlayers
    ind = find(Sites(:,3)==Z(i));
    IND = [IND ind(:)'];
    X = fix((Sites(ind,1)+0.005)*100);
    Y = fix((Sites(ind,2)+0.005)*100);

    for ii=1:length(ind)
        M(X(ii)+9,Y(ii)+8+20*(i-1)) = Qbetadesync(ind(ii));
    end;
end;
imagesc(M);
colorbar
title('Beta desync');

filename_map_beta = sprintf('MapsBeta_%dto%dvs%dto%d_%s.bin',F_bnd1_0,F_bnd1_1,F_bnd2_0,F_bnd2_1, DataFileName{1}(10)); 
fname = sprintf('C:\\LocalData\\Biomag2010Competition\\dataset2preprocessed\\Results\\%s',filename_map_beta);
fid = fopen(fname,'wb');
temp = zeros(1,Nsites*3);
temp(1:3:end) = Qbetadesync;
fwrite(fid,temp,'double');
fclose(fid);

CC = eye(Nsites,Nsites);
%Compute cross terms
P1Pre = zeros(3,Nch);
P2Pre = zeros(3,Nch);
P1Post = zeros(3,Nch);
P2Post = zeros(3,Nch);
for r=1:Nsites
    P1pre =  W1_pre{r}*C1_pre;
    P1post =  W1_post{r}*C1_post;
    P2pre =  W2_pre{r}*C2_pre;
    P2post =  W2_post{r}*C2_post;
    
    for i = r+1:Nsites
        A = P1pre* W1_pre{i}';
        [~, s, ~] = svd(A);
        S1_pre = s(1,1);

        A = P1post* W1_post{i}';
        [~ ,s, ~] = svd(A);
        S1_post = s(1,1);

        A =P2pre*W2_pre{i}';
        [~, s,~] = svd(A);
        S2_pre = s(1,1);

        A = P2post* W2_post{i}';
        [~,s, ~] = svd(A);
        S2_post = s(1,1);
        
        CC(r,i) = (S2_post/S1_post)/(S2_pre/S1_pre);
        CC(i,r) = CC(r,i);
    end
    r
end;

maxx = max(CC,[],2);
IntInd = find(maxx>0.8*max(maxx));
NIntInd = length(IntInd);

% do statistical testing here
tic
R1UpdateBnd1 = [R1UpdatePreBnd1 R1UpdatePostBnd1];
R1UpdateBnd2 = [R1UpdatePreBnd2 R1UpdatePostBnd2];
Npre = size(R1UpdatePreBnd1,2);
Npost = size(R1UpdatePostBnd1,2);
N = Npre+Npost;

for mc = 1:Nmc
    
    C2_post = zeros(Nch);
    C1_post = zeros(Nch);
    C2_pre  = zeros(Nch);
    C1_pre  = zeros(Nch);
    
    if(PermutePreVsPost>0)
        %post
        RandInd = randperm(N);
        indPre = RandInd(1:Npre);

        C1_pre = (R1UpdateBnd1(:,indPre)*R1UpdateBnd1(:,indPre)'+R1UpdateBnd2(:,indPre)*R1UpdateBnd2(:,indPre)')/length(indPre);
        C2_pre = (R1UpdateBnd1(:,indPre)+R1UpdateBnd2(:,indPre))*(R1UpdateBnd1(:,indPre)+R1UpdateBnd2(:,indPre))'/length(indPre);

        indPost = RandInd(Npre+1:end);

        C1_post = (R1UpdateBnd1(:,indPost)*R1UpdateBnd1(:,indPost)'+R1UpdateBnd2(:,indPost)*R1UpdateBnd2(:,indPost)')/length(indPost);
        C2_post = (R1UpdateBnd1(:,indPost)+R1UpdateBnd2(:,indPost))*(R1UpdateBnd1(:,indPost)+R1UpdateBnd2(:,indPost))'/length(indPost);
    else
         for d = 1:size(R1UpdatePostBnd1,2)

            XfftMeanBnd1 = R1UpdatePostBnd1(:,d);
            XfftMeanBnd2 = R1UpdatePostBnd2(:,d);

            %permute the way we compute the two covariances
            if(randn>0)
                C2_post = C2_post + (XfftMeanBnd1+XfftMeanBnd2)*(XfftMeanBnd1+XfftMeanBnd2)';
                C1_post = C1_post + XfftMeanBnd1*XfftMeanBnd1'+XfftMeanBnd2*XfftMeanBnd2';
            else
                C2_post = C2_post + XfftMeanBnd1*XfftMeanBnd1'+XfftMeanBnd2*XfftMeanBnd2';
                C1_post = C1_post + (XfftMeanBnd1+XfftMeanBnd2)*(XfftMeanBnd1+XfftMeanBnd2)';
            end;

        end;

        C1_post = C1_post/size(R1UpdatePostBnd1,2);
        C2_post = C2_post/size(R1UpdatePostBnd1,2);

        %pre
        C1_pre = (R1UpdatePreBnd1*R1UpdatePreBnd1'+R1UpdatePreBnd2*R1UpdatePreBnd2')/size(R1UpdatePreBnd1,2);
        C2_pre = (R1UpdatePreBnd1+R1UpdatePreBnd2)*(R1UpdatePreBnd1+R1UpdatePreBnd2)'/size(R1UpdatePreBnd2,2);

    end;
    Nsites   = size(G,2)/3;
    Nch = size(G,1);

    % covariance data are ready, 
    % build beamformers
    alpha = 0.01; % 0.01 is default for tihinv
    iC1_pre =  tihinv(C1_pre);
    iC2_pre =  tihinv(C2_pre);
    iC1_post = tihinv(C1_post);
    iC2_post = tihinv(C2_post);
    
    W1_pre   = cell(Nsites);
    W2_pre   = cell(Nsites);
    W1_post  = cell(Nsites);
    W2_post  = cell(Nsites);
    
    for i = 1:Nsites
        W1_pre{i} = zeros(3,Nch);
        W2_pre{i} = zeros(3,Nch);
        W1_post{i} = zeros(3,Nch);
        W2_post{i} = zeros(3,Nch);
    end

    disp('computing the inverse operator...');
    range = 1:3;
    for i = 1:Nsites
        denum = tihinv(G(:,range)'*iC1_pre*G(:,range));
        num = G(:,range)'*iC1_pre;
        W1_pre{i} = denum*num;

        denum = tihinv(G(:,range)'*iC1_post*G(:,range));
        num = G(:,range)'*iC1_post;
        W1_post{i} = denum*num;

        denum = tihinv(G(:,range)'*iC2_pre*G(:,range));
        num = G(:,range)'*iC2_pre;
        W2_pre{i} = denum*num;

        denum = tihinv(G(:,range)'*iC2_post*G(:,range));
        num = G(:,range)'*iC2_post;
        W2_post{i} = denum*num;
        range= range+3;
    end;
  
    CCmc = eye(NIntInd,Nsites);
    %Compute cross terms
    P1Pre = zeros(3,Nch);
    P2Pre = zeros(3,Nch);
    P1Post = zeros(3,Nch);
    P2Post = zeros(3,Nch);
    
    disp('computing cross-freq coupling...');
    
    for rr=1:NIntInd
        r = IntInd(rr);
        P1pre =  W1_pre{r}*C1_pre;
        P1post =  W1_post{r}*C1_post;
        P2pre =  W2_pre{r}*C2_pre;
        P2post =  W2_post{r}*C2_post;

        for i = 1:Nsites
            A = P1pre* W1_pre{i}';
            [~, s, ~] = svd(A);
            S1_pre = s(1,1);

            A = P1post* W1_post{i}';
            [~ ,s, ~] = svd(A);
            S1_post = s(1,1);

            A =P2pre*W2_pre{i}';
            [~, s,~] = svd(A);
            S2_pre = s(1,1);

            A = P2post* W2_post{i}';
            [~,s, ~] = svd(A);
            S2_post = s(1,1);

            CCmc(rr,i) = (S2_post/S1_post)/(S2_pre/S1_pre);
        end
        DistMax(rr,mc) = max(CCmc(rr,:));
    end;
    mc
end;

figure

for i=1:NIntInd
    mxActual = max(CC(IntInd(i),:));
    pval(i) = length(find(DistMax(i,:) > mxActual))/Nmc;
    [h,a] = hist(DistMax(i,:),50);
    h = h/max(h);
    subplot(NIntInd,1,i);
    bar(a,h);
    hold on;
    plot([mxActual mxActual],[0 1],'r');
end;
toc

IntInd(find(pval>0.05))=[];
% display significant across-sites interactions

figure
colorbar
K = 150;
Ninterest = length(IntInd);

for k = 1:Ninterest
    M = 0.5*ones(20,20*Nlayers);    
    Xref = fix((Sites(IntInd(k),1))*100);
    Yref = fix((Sites(IntInd(k),2))*100);
    Zref = Sites(IntInd(k),3);
    RefLayerInd = find(Zref==Z);
    IND = [];
    for i=1:Nlayers
        ind = find(Sites(:,3)==Z(i));
        IND = [IND ind(:)'];
        X = fix((Sites(ind,1)+0.005)*100);
        Y = fix((Sites(ind,2)+0.005)*100);
        
        for ii=1:length(ind)
            M(X(ii)+9,Y(ii)+8+20*(i-1)) = CC(IntInd(k),ind(ii));
        end;
        if(i==RefLayerInd)
            tmp = max(CC(IntInd(k),:));
            M(Xref+9,Yref+8+20*(i-1)) = tmp;
            M(Xref+9-1,Yref+8+20*(i-1)) = tmp;
            M(Xref+9-1,Yref+8+20*(i-1)+1) = tmp;
            M(Xref+9,Yref+8+20*(i-1)+1) = tmp;
        end;
    end;
    subplot(Ninterest,1,k);
    imagesc(M);
    colorbar
end;

if(PermutePreVsPost)
    filename_map = sprintf('Maps_%dto%dvs%dto%d_%s_PreVsPost_%d',F_bnd1_0,F_bnd1_1,F_bnd2_0,F_bnd2_1, DataFileName{1}(10),PermutePreVsPost); 
else
    filename_map = sprintf('Maps_%dto%dvs%dto%d_%s_FullVsCodeprived_%d',F_bnd1_0,F_bnd1_1,F_bnd2_0,F_bnd2_1, DataFileName{1}(10),PermutePreVsPost); 
end;


%save significant ones to visualize in EMSE
for i=1:Ninterest
    fname = sprintf('C:\\LocalData\\Biomag2010Competition\\dataset2preprocessed\\Results\\%s_N_%d.bin',filename_map,i);
    fid = fopen(fname,'wb');
    temp = zeros(1,Nsites*3);
    V = CC(IntInd(i),:);
    V(IntInd(i)) = max(V);
    temp(1:3:end) = V;
    fwrite(fid,temp,'double');
    fclose(fid);
    i
end;

