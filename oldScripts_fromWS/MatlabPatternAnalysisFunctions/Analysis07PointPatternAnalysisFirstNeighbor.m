function Analysis03PointPatternAnalysisFirstNeighbor(path,filename,k,x,y,z,S,d123_all,d123_1)
%Point pattern analysis Type 2 distance to Type 3 nearest cells in Type 1+2 (first neighbor)
%Modified September 27, 2017 

colormap(jet);
CC=colormap;



S=S';

CellDiameter=mean(d123_1);%check this definition !


index2=find(S==2);
index3=find(S==3);
for i=1:length(index2);
    dtemp=sort(d123_all(index2(i),index3));
    dn(i)=dtemp(1);%nearest Type 3 to each Type 2 cell
end
dn=dn/CellDiameter;
dndata=dn;
r=0:0.1:14;
for i=1:size(r,2)
    G(i)=size(dn(dn<=r(i)),2)/size(dn,2);
end
figure
plot(r,G,'-r','linewidth',4);hold on
ylabel('Cumulative cell frequency');
xlabel('distance to nearest Type 3 neighbor to Type 2 cells (in cell diameter)');



if 2*ceil(r(min(find(G>0.9))))<14
    axis([0 2*ceil(r(min(find(G>0.9)))) 0 1])
else
    axis([0 14 0 1])
end

%         xin=1:1:10;
ii=0;
D=[];
NumPermut=1000;
d1min=0;%in µm
for j=1:NumPermut
    j

    index12=find(S==1|S==2);
    S(index12)=S(index12(randperm(length(index12))));

    if d1min>0
        [dtemp d1 d2]=NNanalysis(x(S==2),y(S==2),z(S==2));
        while min(d1)<d1min
            S(index12)=S(index12(randperm(length(index12))));
            [dtemp d1 d2]=NNanalysis(x(S==2),y(S==2),z(S==2));
        end
        clear dn dtemp
    end
    clear dn dtemp
    index2=find(S==2);
    index3=find(S==3);

    for i=1:length(index2);
        dtemp=sort(d123_all(index2(i),index3));
        dn(i)=dtemp(1);%nearest Type 3 to each Type 2 cell
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




for i=1:size(r,2)
    Grandmean(i)=mean(Grand(:,i));
    Grandstd(i)=std(Grand(:,i));
    Grand5(i)=prctile(Grand(:,i),5);
    Grand95(i)=prctile(Grand(:,i),95);
    Grand1(i)=prctile(Grand(:,i),1);
    Grand99(i)=prctile(Grand(:,i),99);
end
plot(r,Grandmean,'-k','linewidth',4);hold on
%     plot(r,Grandmean-Grandstd,'-g');hold on
%     plot(r,Grandmean+Grandstd,'-g');hold on
plot(r,Grand5,'-','linewidth',3,'color',[0.6 0.6 0.6]);hold on
plot(r,Grand95,'-','linewidth',3,'color',[0.6 0.6 0.6]);hold on
plot(r,Grand1,'-','linewidth',1,'color',[0.6 0.6 0.6]);hold on
plot(r,Grand99,'-','linewidth',1,'color',[0.6 0.6 0.6]);hold on


plot(r,G,'-r','linewidth',4);hold on
text(0.5,0.95,['95% and 99% intervals' ])
text(0.5,0.9,[num2str(NumPermut),' random perm.'])
%         text(0.5,0.85,[num2str(d1min),'µm min dist.' ])
title(['Type II (random permutations in I and II) CellDiameter=',num2str(CellDiameter,2),'µm'],'FontSize',10)
saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis07Fig1'], 'tiffn');
clear dn
clear G









