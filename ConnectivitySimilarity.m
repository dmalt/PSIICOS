function CS = ConnectivitySimilarity(Pairs1,Pairs2,ChLoc, NOrder,MaxDist)

if(nargin==3)
    NOrder = 4;
    MaxDist = 0.05;
end;

for i=1:size(ChLoc,2)
    DeltaVec  = bsxfun(@minus, ChLoc,ChLoc(:,i));
    Delta = sqrt(sum(DeltaVec.^2,1));
    [DeltaSrt,Key_srt] = sort(Delta);
    Neighs = Key_srt(2:(NOrder+1));
    Neighs(DeltaSrt(2:(NOrder+1))>MaxDist) = [];
    ElNeighbs{i} = Neighs;
end;
M1 = zeros(204,204);
for i=1:length(Pairs1)
    M1(Pairs1(i,1),Pairs1(i,2)) = 1;
    M1(Pairs1(i,2),Pairs1(i,1)) = 1;
    for n1 = 1:length(ElNeighbs{Pairs1(i,1)})
        for n2 = 1:length(ElNeighbs{Pairs1(i,2)})
            M1(ElNeighbs{Pairs1(i,1)}(n1),ElNeighbs{Pairs1(i,2)}(n2)) = 1;
        end;
    end;
end;

M2 = zeros(204,204);
for i=1:length(Pairs2)
    M2(Pairs2(i,1),Pairs2(i,2)) = 1;
    M2(Pairs2(i,2),Pairs2(i,1)) = 1;    
    for n1 =length(ElNeighbs{Pairs2(i,1)})
        for n2 =length( ElNeighbs{Pairs2(i,2)})
            M2(ElNeighbs{Pairs2(i,1)}(n1),ElNeighbs{Pairs2(i,2)}(n2)) = 1;
        end;
    end;
end;

CS = sum(sum(M1.*M2))/(0.5*(length(Pairs1)+length(Pairs2)));

