function fullResults = analysis05PPAFirstNeighbor(path,name,k,x,y,z,S,d123_all,d123_1)
%Point pattern analysis Type 3 in Type 1+2+3 (first neighbor)
r=0:0.1:14; % bin size for the ecdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Type 1+2+3

popSource = 3;
popTarget = 3;
popPermut = [1 2 3];

S=S';
CellDiameter=mean(d123_1);

index=find(S==3);
for i=1:length(S(S==3))
    dtemp=sort(d123_all(index(i),index));
    dn(i)=dtemp(2);
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
d1min=0;%in um
for j=1:NumPermut % for each permutation run
    S=S(randperm(length(S)));
    if d1min>0
        [dtempSimu d1 d2]=NNanalysis(x(S==popTarget),y(S==popTarget),z(S==popTarget));
        while min(d1)<d1min
            S=S(randperm(length(S)));
            [dtempSimu d1 d2]=NNanalysis(x(S==popTarget),y(S==popTarget),z(S==popTarget));
        end
        clear dtempSimu
    end
    %     clear dtemp
    index=find(S==popTarget);
    
    for i=1:numel(index)
        dtempSimu=sort(d123_all(index(i),index));
        dnSimu(i)=dtempSimu(2);
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
figSavePath = [path,name,'Case',num2str(k),'_Analysis05Fig1'];

GrandAll = displaySimuPPQ(r, Grand, G, NumPermut, figTitle, figSavePath);

save([path,name,'Case',num2str(k),'_Analysis05NN'],'dn','G','r','Grand','GrandAll');


fullResults = {};
fullResults.dn = dn;
fullResults.G = G;
fullResults.GrandAll = GrandAll;

end









