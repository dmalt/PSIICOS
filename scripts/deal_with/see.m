Band = [18,21];
Fsamp = 500;
range = 201:250;
ConData;

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
