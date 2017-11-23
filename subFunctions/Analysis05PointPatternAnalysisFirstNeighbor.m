function Analysis05PointPatternAnalysisFirstNeighbor(path,name,k,x,y,z,S,d123_all,d123_1)
%Point pattern analysis Type 3 in Type 1+2+3 (first neighbor)

colormap(jet);
CC=colormap;
% MaxDistance=14;%in CellDiameter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type 1+2+3

%     Snew=S(S==1|S==2);
%     S=Snew';
S=S';
CellDiameter=mean(d123_1);

index=find(S==3);
for i=1:length(S(S==3))
    dtemp=sort(d123_all(index(i),index));
    dn(i)=dtemp(2);
end
dn=dn/CellDiameter;
dndata=dn;
r=0:0.1:14;
for i=1:size(r,2)
    G(i)=size(dn(dn<=r(i)),2)/size(dn,2);
end

if 2*ceil(r(min(find(G>0.9))))<14
    axis([0 2*ceil(r(min(find(G>0.9)))) 0 1])
else
    axis([0 14 0 1])
end

%         xin=1:1:10;
ii=0;
D=[];
NumPermut=1000;
d1min=0;%in �m
for j=1:NumPermut
    S=S(randperm(length(S)));
    if d1min>0
        [dtemp d1 d2]=NNanalysis(x(S==3),y(S==3),z(S==3));
        while min(d1)<d1min
            S=S(randperm(length(S)));
            [dtemp d1 d2]=NNanalysis(x(S==3),y(S==3),z(S==3));
        end
        clear dn dtemp
    end
    clear dn dtemp
    index=find(S==3);
    
    for i=1:length(S(S==3))
        dtemp=sort(d123_all(index(i),index));
        dn(i)=dtemp(2);
    end
    dn=dn/CellDiameter;
    %             r=0:0.1:10;
    %             [n(j,:),xout] = hist(dn,xin);
    for i=1:size(r,2)
        Grand(j,i)=size(dn(dn<=r(i)),2)/size(dn,2);
    end
end

%         save([name,'ALLrandperm',num2str(NumPermut),'d1min',num2str(d1min),'.mat'],'Grand');
%         load([folder,name,'ALLrandperm',num2str(NumPermut),'d1min',num2str(d1min),'.mat'])

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









