function Analysis01DistributionDisplay(path,filename,k,x,y,z,S,d1meanType123,d1meanType12,d1meanType13,ColorRange)
%Modified April 21, 2017 

colormap(jet);
CC=colormap;
    
    figure
    set(gcf,'Position',[20 20 1000 800]);
    plot(x(S==2),y(S==2),'o','LineWidth',2,'MarkerSize',2,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);hold on;
    plot(x(S==1),y(S==1),'o','LineWidth',1,'MarkerSize',1,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);hold on;
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold');
    title(['N_{I+II}=',num2str(length(S(S==1|S==2))),' N_{II}=',num2str(length(S(S==2))),' (',num2str(100*length(S(S==2))/length(S(S==1|S==2)),2),'%)',' <d1Type1&2>=',num2str(d1meanType12,2),'µm'])
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis01Fig1'], 'tiffn');
    
    
    figure
    set(gcf,'Position',[20 20 1000 800]);
   plot(x(S==3),y(S==3),'o','LineWidth',2,'MarkerSize',2,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);hold on;
    plot(x(S==1),y(S==1),'o','LineWidth',1,'MarkerSize',1,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);hold on;
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold');
    title(['N_{I+III}=',num2str(length(S(S==1|S==3))),' N_{III}=',num2str(length(S(S==3))),' (',num2str(100*length(S(S==3))/length(S(S==1|S==3)),2),'%)',' <d1Type1&3>=',num2str(d1meanType13,2),'µm'])
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis01Fig2'], 'tiffn');
%    
%     figure
%     set(gcf,'Position',[20 20 1000 800]);
%     % plot3(x(S==1),y(S==1),z(S==1),'o','LineWidth',2,'MarkerSize',3,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);hold on;
%     % plot3(x(S==0),y(S==0),z(S==0),'o','LineWidth',1,'MarkerSize',2,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);hold on;
%     plot(x(S==3),y(S==3),'o','LineWidth',2,'MarkerSize',2,'MarkerFaceColor',[0 0 1],'MarkerEdgeColor',[0 0 1]);hold on;
%     plot(x(S==2),y(S==2),'o','LineWidth',2,'MarkerSize',2,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);hold on;
%     plot(x(S==1),y(S==1),'o','LineWidth',1,'MarkerSize',1,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);hold on;
%     axis equal;
%     axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
%     %  axis vis3d;
%     set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold');
%     % view(45,-50)
%     % if k==5
%     %     view(-61,-86)
%     % elseif k==4
%     %     view(43,90)
%     % end
%     % zoom(1.45)
%     
%     title([name,' N_{I+II+III}=',num2str(size(x,1)),' N_{II}=',num2str(length(S(S==2))),' (',num2str(100*length(S(S==2))/size(x,1),2),'%)',' <dn>=',num2str(dmean,2),'µm'])
%     saveas(gcf,['Analysis04',name,'Fig3'], 'tiffn');
    
    
    [d d1 d2]=NNanalysis(x(S==2),y(S==2),z(S==2));
    d4=(d1-ColorRange(1))/(ColorRange(2)-ColorRange(1));    
    figure
    set(gcf,'Position',[20 20 1000 800]);
  
    index=find(S==2);
    for i=1:length(index)

        ii=1+floor((1-1*d4(i))*63);
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        plot(x(index(i)),y(index(i)),'o','LineWidth',2,'MarkerSize',4+ii/64*2,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold');
    title(['Type II nearest neighbor distance range=[',num2str(ColorRange(1)),'µm ',num2str(ColorRange(2)),'µm] mean=',num2str(mean(d1),2),'µm'],'FontSize',10)
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis01Fig3'], 'tiffn');
    
    
    
    [d d1 d2]=NNanalysis(x(S==3),y(S==3),z(S==3));
    d4=(d1-ColorRange(1))/(ColorRange(2)-ColorRange(1));
    
    
        figure
    set(gcf,'Position',[20 20 1000 800]);
    index=find(S==3);
    for i=1:length(index)

        ii=1+floor((1-1*d4(i))*63);
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        plot(x(index(i)),y(index(i)),'o','LineWidth',2,'MarkerSize',4+ii/64*2,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    set(gca,'LineWidth',1.5,'FontSize',12,'FontWeight','bold');
    title(['Type III nearest neighbor distance range=[',num2str(ColorRange(1)),'µm ',num2str(ColorRange(2)),'µm] mean=',num2str(mean(d1),2),'µm'],'FontSize',10)
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis01Fig4'], 'tiffn');
      

