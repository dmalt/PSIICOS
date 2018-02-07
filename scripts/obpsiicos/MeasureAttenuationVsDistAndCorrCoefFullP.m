function [NormAre, NormArep, NormAim, NormAimp,NormVC1, NormVC1p,NormVC2, NormVC2p,RFirstSecond, CFirstSecond] = MeasureAttenuationVsDistAndCorrCoefFullP(Pr,G,R,N,Nmc)

Nsrc = size(G,2)/2;

if(N>Nsrc)
    disp('Number of pairs can not be greater than half the number of the forward matrix columns.');
    return 
end;

RandSrcInds = randperm(Nsrc*2);
FirstSrcInds = RandSrcInds(1:N);
SecondSrcInds = RandSrcInds(N+1:N+N);
SourcePairs = [FirstSrcInds(:) SecondSrcInds(:)];

Nch  = size(G,1);
Are = zeros(Nch*Nch,size(SourcePairs,1));
Aim = zeros(Nch*Nch,size(SourcePairs,1));
Avc1 = zeros(Nch*Nch,size(SourcePairs,1));
Avc2 = zeros(Nch*Nch,size(SourcePairs,1));

Np = size(SourcePairs,1);
NormAre = zeros(1,N*Nmc);
NormAim = zeros(1,N*Nmc);
NormArep = zeros(1,N*Nmc);
NormAimp = zeros(1,N*Nmc);
NormVC1 = zeros(1,N*Nmc);
NormVC1p = zeros(1,N*Nmc);
NormVC2 = zeros(1,N*Nmc);
NormVC2p = zeros(1,N*Nmc);
RFirstSecond = zeros(1,N*Nmc);
CFirstSecond = zeros(1,N*Nmc);

P = 0;
for mc = 1:Nmc
    
    RandSrcInds   = randperm(Nsrc*2);
    FirstSrcInds  = RandSrcInds(1:N);
    SecondSrcInds = RandSrcInds(N+1:N+N);
    SourcePairs   = [FirstSrcInds(:) SecondSrcInds(:)];

    for p=1:Np
        i1 = SourcePairs(p,1);
        i2 = SourcePairs(p,2);
        tmp1 = G(:,i1)*G(:,i2)';
        tmp2 = tmp1+tmp1';
        Are(:,p) = tmp2(:);
        tmp2 = tmp1-tmp1';
        Aim(:,p) = tmp2(:);
        tmp1 = G(:,i1)*G(:,i1)';
        tmp2 = G(:,i2)*G(:,i2)';
        Avc1(:,p) = tmp1(:); 
        Avc2(:,p) = tmp2(:);
        RFirstSecond(P+p) = norm(R(i1,:)-R(i2,:));
        CFirstSecond(P+p) = (G(:,i1)'*G(:,i2))/(norm(G(:,i1))*norm(G(:,i2)));
    end;

    Arep  = Pr*Are; %Are-Upwr*(Upwr'*Are);
    Aimp  = Pr*Aim; %-Upwr*(Upwr'*Aim);
    Avc1p  = Pr*Avc1;%-Upwr*(Upwr'*Avc1);
    Avc2p  = Pr*Avc2;%-Upwr*(Upwr'*Avc2);
    
    
    NormAre(P+1:P+N) = sqrt(sum(Are.*Are,1));
    NormAim(P+1:P+N) = sqrt(sum(Aim.*Aim,1));
    NormVC1(P+1:P+N) = sqrt(sum(Avc1.*Avc1,1));
    NormVC2(P+1:P+N) = sqrt(sum(Avc2.*Avc2,1));
  
    NormArep(P+1:P+N) = sqrt(sum(Arep.*Arep,1));
    NormAimp(P+1:P+N) = sqrt(sum(Aimp.*Aimp,1));
    NormVC1p(P+1:P+N) = sqrt(sum(Avc1p.*Avc1p,1));
    NormVC2p(P+1:P+N) = sqrt(sum(Avc2p.*Avc2p,1));
        
    P = P+N;
end; %mc


