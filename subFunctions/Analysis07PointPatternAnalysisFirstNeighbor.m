function Analysis07PointPatternAnalysisFirstNeighbor(path,name,k,x,y,z,S,d123_all,d123_1)
%Point pattern analysis Type 2 distance to Type 3 nearest cells in Type 1+2 (first neighbor)
% Could use a list of all diameters instead of sending the specific one...
r=0:0.1:14; % bin size for the ecdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type 1+2+3

popTarget = 2; % pop of interest for the distribution
popSource = 3; % pop to which the target pop distribution is compared
popPermut = [1 2]; % pop used for the random permutation

S=S';
CellDiameter=mean(d123_1);


% Find nearest neighbour from popTarget to popSource
indexPopTarget = find(S==popTarget);
indexPopSource = find(S==popSource);
for i=1:length(indexPopTarget)
    dtemp=sort(d123_all(indexPopTarget(i),indexPopSource));
    dn(i)=dtemp(1);%nearest Type 3 to each Type 2 cell
end

% Normalize min distances by cell Diameter (using 123 or 12 only?)
dn=dn/CellDiameter;

% calculate the cdf by hand... consider using ecdf for proper method
for i=1:size(r,2)
    G(i)=size(dn(dn<=r(i)),2)/size(dn,2);
end


% Effectuate the simulation of the random permutations
ii=0;
D=[];
NumPermut=1000;
d1min=0;%in �m
for j=1:NumPermut % for each permutation run
    index12=find(S==1|S==2);
    S(index12)=S(index12(randperm(length(index12))));
    if d1min>0
        [dtempSimu d1 d2]=NNanalysis(x(S==popTarget),y(S==popTarget),z(S==popTarget));
        while min(d1)<d1min
            S(index12)=S(index12(randperm(length(index12))));
            [dtempSimu d1 d2]=NNanalysis(x(S==popTarget),y(S==popTarget),z(S==popTarget));
        end
        clear dtempSimu
    end
    %     clear dtemp
    indexPopTargetSimu=find(S==popTarget);
    indexPopSourceSimu=find(S==popSource);
    
    for i=1:length(indexPopTargetSimu)
        dtempSimu=sort(d123_all(indexPopTargetSimu(i),indexPopSourceSimu));
        dnSimu(i)=dtempSimu(1);%nearest Type 3 to each Type 2 cell
    end
    dnSimu=dnSimu/CellDiameter;
    %             r=0:0.1:10;
    %             [n(j,:),xout] = hist(dn,xin);
    for i=1:size(r,2)
        Grand(j,i)=size(dnSimu(dnSimu<=r(i)),2)/size(dnSimu,2);
    end
end

% Display figure
figTitle = ['Type II (random permutations in I and II) CellDiameter=',num2str(CellDiameter,2),'{\mu}m'];
figSavePath = [path,name,'Case',num2str(k),'_Analysis07Fig1'];

GrandAll = displaySimuPPQ(r, Grand, G, NumPermut, figTitle, figSavePath);

Grandmean = GrandAll.mean;
Grandstd = GrandAll.std;
Grand5 = GrandAll.iqr5;
Grand95 = GrandAll.iqr95;
Grand1 = GrandAll.iqr1;
Grand99 = GrandAll.iqr99;

save([path,name,'Case',num2str(k),'_Analysis07NN'],'dn','G','r','Grand','Grandmean','Grandstd',...
    'Grand5','Grand95','Grand1','Grand99');

clear dn
clear G









