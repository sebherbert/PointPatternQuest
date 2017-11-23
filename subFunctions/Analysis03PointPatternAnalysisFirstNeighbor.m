function Analysis03PointPatternAnalysisFirstNeighbor(path,name,k,x,y,z,S,d12_all,d12_1)
%Point pattern analysis Type 2 in Type 1+2 (first neighbor)

colormap(jet);
CC=colormap;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type 1+2 only

Snew=S(S==1|S==2);
S=Snew';
CellDiameter=mean(d12_1);

index=find(S==2);
for i=1:length(S(S==2))
    dtemp=sort(d12_all(index(i),index));
    dn(i)=dtemp(2);
end
dn=dn/CellDiameter;
r=0:0.1:14;
for i=1:size(r,2)
    G(i)=size(dn(dn<=r(i)),2)/size(dn,2);
end

if 2*ceil(r(min(find(G>0.9))))<14
    axis([0 2*ceil(r(min(find(G>0.9)))) 0 1]);
else
    axis([0 14 0 1]);
end

%         xin=1:1:10;
ii=0;
D=[];
NumPermut=1000;
d1min=0;%in �m
for j=1:NumPermut
    S=S(randperm(length(S)));
    if d1min>0
        [dtemp d1 d2]=NNanalysis(x(S==2),y(S==2),z(S==2));
        while min(d1)<d1min
            S=S(randperm(length(S)));
            [dtemp d1 d2]=NNanalysis(x(S==2),y(S==2),z(S==2));
        end
        clear dn dtemp
    end
    clear dn dtemp
    index=find(S==2);
    
    for i=1:length(S(S==2))
        dtemp=sort(d12_all(index(i),index));
        dn(i)=dtemp(2);
    end
    dn=dn/CellDiameter;
    %             r=0:0.1:10;
    %             [n(j,:),xout] = hist(dn,xin);
    for i=1:size(r,2)
        Grand(j,i)=size(dn(dn<=r(i)),2)/size(dn,2);
    end
end

% Display figure
figTitle(['Type II (random permutations in I and II) CellDiameter=',num2str(CellDiameter,2),'{\mu}m'],'FontSize',10);
figSavePath([path,name,'Case',num2str(k),'_Analysis07Fig1']);

GrandAll = displaySimuPPQ(r, Grand, G, NumPermut, figTitle, figSavePath);

Grandmean = GrandAll.mean;
Grandstd = GrandAll.std;
Grand5 = GrandAll.iqr5;
Grand95 = GrandAll.iqr95;
Grand1 = GrandAll.iqr1;
Grand99 = GrandAll.iqr99;

save([path,name,'Case',num2str(k),'_Analysis05NN'],'dn','G','r','Grand','Grandmean','Grandstd',...
    'Grand5','Grand95','Grand1','Grand99');

clear dn
clear G









