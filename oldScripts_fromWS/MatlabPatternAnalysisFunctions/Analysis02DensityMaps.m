function Analysis021DensityMaps(path,filename,k,x,y,z,S,d12_all,d123_all,d12_1,d123_1,LocalRadius,ColorRange)
%Modified July 12, 2017 
colormap(jet);
CC=colormap;
    
ColorRange(5,:)=ColorRange(5,:)/100;
ColorRange(6,:)=ColorRange(6,:)/100;
ColorRange(7,:)=ColorRange(7,:)/100;
ColorRange(8,:)=ColorRange(8,:)/100;
ColorRange(9,:)=ColorRange(9,:)/100;

   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fig1: Local Density Type 1+2+3
        
    densityType2mean=length(S(S==2))/length(x);
    densityType3mean=length(S(S==3))/length(x);
    densityType2and3mean=length(S(S==3|S==2))/length(x);

     
    for i=1:size(x);    
        index=find(d123_all(i,:)<LocalRadius);%local density in a LocalRadius µm area
        density(i)=mean(d123_1(index));
        
        N1=sum(S(index)==1);
        N2=sum(S(index)==2);
        N3=sum(S(index)==3);
        
        densityType2(i)=N2/(N1+N2+N3);
        densityType3(i)=N3/(N1+N2+N3);
        densityType2and3(i)=(N2+N3)/(N1+N2+N3);
    end
    
    figure
    set(gcf,'Position',[100 100 1000 800]);
    for i=1:size(x)
        ii=1+floor(63-(density(i)-ColorRange(1,1))/(ColorRange(1,2)-ColorRange(1,1))*63);
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        plot(x(i),y(i),'o','LineWidth',1,'MarkerSize',1+ii/64*4,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    title(['Local density I+II+III (local average distance to nearest neighbor) range=[',num2str(ColorRange(1,1)),'µm ',num2str(ColorRange(1,2)),'µm] mean=',num2str(mean(d123_1),2),'µm'],'FontSize',10)
    set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');    
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis02Fig1'], 'tiffn');
    
    xlsfilename=[path,filename(1:end-4),'Case',num2str(k),'_Analysis02'];
    delete([xlsfilename,'.xls'])
    %Write to excel file
    xlspage='Fig1';
    xlswrite(xlsfilename, {'x','y','z','Local density I+II+III'}, xlspage, 'A1');
    xlswrite(xlsfilename,x,xlspage,'A2');
    xlswrite(xlsfilename,y,xlspage,'B2');
    xlswrite(xlsfilename,z,xlspage,'C2');
    xlswrite(xlsfilename,density',xlspage,'D2');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fig5: local ratio II/I+II+III   
    figure
    set(gcf,'Position',[100 100 1000 800]);
    for i=1:size(x)
        ii=1+floor((densityType2(i)-ColorRange(5,1))/(ColorRange(5,2)-ColorRange(5,1))*63);         
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        plot(x(i),y(i),'o','LineWidth',1,'MarkerSize',1+ii/64*4,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    title(['Local ratio II/(I+II+III) density range=[',num2str(ColorRange(5,1)*100),'% ',num2str(ColorRange(5,2)*100),'%] mean=',num2str(densityType2mean*100,2),'%'],'FontSize',10)
    set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis02Fig5'], 'tiffn');
    
    %Write to excel file
    xlspage='Fig5';
    xlswrite(xlsfilename, {'x','y','z','Local ratio II/I+II+III'}, xlspage, 'A1');
    xlswrite(xlsfilename,x,xlspage,'A2');
    xlswrite(xlsfilename,y,xlspage,'B2');
    xlswrite(xlsfilename,z,xlspage,'C2');
    xlswrite(xlsfilename,densityType2',xlspage,'D2');
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fig6: local ratio III/I+II+III
    figure
    set(gcf,'Position',[100 100 1000 800]);
    for i=1:size(x)
        ii=1+floor((densityType3(i)-ColorRange(6,1))/(ColorRange(6,2)-ColorRange(6,1))*63);         
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        plot(x(i),y(i),'o','LineWidth',1,'MarkerSize',1+ii/64*4,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    title(['Local ratio III/(I+II+III) density range=[',num2str(ColorRange(6,1)*100),'% ',num2str(ColorRange(6,2)*100),'%] mean=',num2str(densityType3mean*100,2),'%'],'FontSize',10)
    set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis02Fig6'], 'tiffn');
    
    %Write to excel file
    xlspage='Fig6';
    xlswrite(xlsfilename, {'x','y','z','Local ratio III/I+II+III'}, xlspage, 'A1');
    xlswrite(xlsfilename,x,xlspage,'A2');
    xlswrite(xlsfilename,y,xlspage,'B2');
    xlswrite(xlsfilename,z,xlspage,'C2');
    xlswrite(xlsfilename,densityType3',xlspage,'D2');
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fig7: local ratio II+III/I+II+III    
    figure
    set(gcf,'Position',[100 100 1000 800]);
    for i=1:size(x)
        ii=1+floor((densityType2and3(i)-ColorRange(7,1))/(ColorRange(7,2)-ColorRange(7,1))*63);
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        plot(x(i),y(i),'o','LineWidth',1,'MarkerSize',1+ii/64*4,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    title(['Local ratio (II+III)/(I+II+III) density range=[',num2str(ColorRange(7,1)*100),'% ',num2str(ColorRange(7,2)*100),'%] mean=',num2str(densityType2and3mean*100,2),'%'],'FontSize',10)
    set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis02Fig7'], 'tiffn');
    
    %Write to excel file
    xlspage='Fig7';
    xlswrite(xlsfilename, {'x','y','z','Local ratio II+III/I+II+III'}, xlspage, 'A1');
    xlswrite(xlsfilename,x,xlspage,'A2');
    xlswrite(xlsfilename,y,xlspage,'B2');
    xlswrite(xlsfilename,z,xlspage,'C2');
    xlswrite(xlsfilename,densityType2and3',xlspage,'D2');

    
    clear whos density*
    clear N1 N2 N3
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fig2: Local Density Type 1+2
      
    densityType2mean=length(S(S==2))/length(x(S==1|S==2));
    Snew=S(S==1|S==2);
    
    for i=1:length(x(S==1|S==2));    
        index=find(d12_all(i,:)<LocalRadius);%local density in a LocalRadius µm area
        density(i)=mean(d12_1(index));

        N1=sum(Snew(index)==1);
        N2=sum(Snew(index)==2);
        densityType2(i)=N2/(N1+N2);
    end
    
    
    
    figure
    set(gcf,'Position',[100 100 1000 800]);
    for i=1:length(x(S==1|S==2))
        ii=1+floor(63-(density(i)-ColorRange(2,1))/(ColorRange(2,2)-ColorRange(2,1))*63);
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        plot(x(i),y(i),'o','LineWidth',1,'MarkerSize',1+ii/64*4,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    title(['Local density I+II (local average distance to nearest neighbor) range=[',num2str(ColorRange(2,1)),'µm ',num2str(ColorRange(2,2)),'µm] mean=',num2str(mean(d12_1),2),'µm'],'FontSize',10)
    set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');    
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis02Fig2'], 'tiffn');
        
    %Write to excel file
    xlspage='Fig2';
    xlswrite(xlsfilename, {'x','y','z','Local Density Type 1+2'}, xlspage, 'A1');
    xlswrite(xlsfilename,x(S==1|S==2),xlspage,'A2');
    xlswrite(xlsfilename,y(S==1|S==2),xlspage,'B2');
    xlswrite(xlsfilename,z(S==1|S==2),xlspage,'C2');
    xlswrite(xlsfilename,density',xlspage,'D2');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fig8: local ratio II/I+II 
    figure
    set(gcf,'Position',[100 100 1000 800]);
    for i=1:length(x(S==1|S==2))
        ii=1+floor((densityType2(i)-ColorRange(8,1))/(ColorRange(8,2)-ColorRange(8,1))*63);         
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        plot(x(i),y(i),'o','LineWidth',1,'MarkerSize',1+ii/64*4,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    title(['Local ratio II/(I+II) density range=[',num2str(ColorRange(8,1)*100),'% ',num2str(ColorRange(8,2)*100),'%] mean=',num2str(densityType2mean*100,2),'%'],'FontSize',10)
    set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis02Fig8'], 'tiffn');
    
          
    %Write to excel file
    xlspage='Fig8';
    xlswrite(xlsfilename, {'x','y','z','Local ratio II/I+II'}, xlspage, 'A1');
    xlswrite(xlsfilename,x(S==1|S==2),xlspage,'A2');
    xlswrite(xlsfilename,y(S==1|S==2),xlspage,'B2');
    xlswrite(xlsfilename,z(S==1|S==2),xlspage,'C2');
    xlswrite(xlsfilename,densityType2',xlspage,'D2');
    

    
    clear whos density*
    clear N1 N2 N3
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FigX: Local Density Type 2+3
      
    densityType2mean=length(S(S==2))/length(x(S==2|S==3));
    Snew=S(S==2|S==3);
    [d23_all,d23_1,d23_2]=NNanalysis(x(S==2|S==3),y(S==2|S==3),z(S==2|S==3));
    
    
    inew=find(S==2|S==3);
    for i=1:length(x(S==2|S==3));
        
        index=find(d23_all(i,:)<LocalRadius);%local density in a LocalRadius µm area
        density(i)=mean(d23_1(index));

        N2=sum(Snew(index)==2);
        N3=sum(Snew(index)==3);
        densityType2(i)=N2/(N2+N3);
    end
    
%     figure
%     set(gcf,'Position',[100 100 1000 800]);
%     for i=1:length(x(S==2|S==3))
%         ii=1+floor(63-(density(i)-ColorRange(X,1))/(ColorRange(X,2)-ColorRange(X,1))*63);
%         if ii<1
%             ii=1;
%         elseif ii>64
%             ii=64;
%         end
%         colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
%         plot(x(i),y(i),'o','LineWidth',1,'MarkerSize',1+ii/64*4,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
%     end
%     axis equal;
%     axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
%     title(['Local density II+II (local average distance to nearest neighbor) range=[',num2str(ColorRange(X,1)),'µm ',num2str(ColorRange(X,2)),'µm] mean=',num2str(mean(d23_1),2),'µm'],'FontSize',10)
%     set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');    
%     saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis02FigX'], 'tiffn');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fig9: local ratio II/II+III 
    figure
    set(gcf,'Position',[100 100 1000 800]);
%     plot(x(S==2),y(S==2),'r+');hold on;
%     plot(x(S==3),y(S==3),'b+');
    for i=1:length(x(S==2|S==3))
        ii=1+floor((densityType2(i)-ColorRange(9,1))/(ColorRange(9,2)-ColorRange(9,1))*63);         
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        
        plot(x(inew(i)),y(inew(i)),'o','LineWidth',1,'MarkerSize',1+ii/64*4,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    title(['Local ratio II/(II+III) density range=[',num2str(ColorRange(9,1)*100),'% ',num2str(ColorRange(9,2)*100),'%] mean=',num2str(densityType2mean*100,2),'%'],'FontSize',10)
    set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis02Fig9'], 'tiffn');
              
    %Write to excel file
    xlspage='Fig9';
    xlswrite(xlsfilename, {'x','y','z','Local ratio II/II+III'}, xlspage, 'A1');
    xlswrite(xlsfilename,x(inew),xlspage,'A2');
    xlswrite(xlsfilename,y(inew),xlspage,'B2');
    xlswrite(xlsfilename,z(inew),xlspage,'C2');
    xlswrite(xlsfilename,densityType2',xlspage,'D2');
    
    
   clear whos density*
    clear N1 N2 N3  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fig3: Local Density Type 2
      
    [d2_all d2_1 d2_2]=NNanalysis(x(S==2),y(S==2),z(S==2));
    
    
    inew=find(S==2);
    for i=1:length(x(S==2));       
        index=find(d2_all(i,:)<LocalRadius);%local density in a LocalRadius µm area
        density(i)=mean(d2_1(index));
    end
    
    figure
    set(gcf,'Position',[100 100 1000 800]);
    for i=1:length(x(S==2))
        ii=1+floor(63-(density(i)-ColorRange(3,1))/(ColorRange(3,2)-ColorRange(3,1))*63);
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        plot(x(inew(i)),y(inew(i)),'o','LineWidth',1,'MarkerSize',1+ii/64*4,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    title(['Local density II (local average distance to nearest neighbor) range=[',num2str(ColorRange(3,1)),'µm ',num2str(ColorRange(3,2)),'µm] mean=',num2str(mean(d2_1),2),'µm'],'FontSize',10)
    set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');    
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis02Fig3'], 'tiffn');
    
    %Write to excel file
    xlspage='Fig3';
    xlswrite(xlsfilename, {'x','y','z','Local Density Type II'}, xlspage, 'A1');
    xlswrite(xlsfilename,x(inew),xlspage,'A2');
    xlswrite(xlsfilename,y(inew),xlspage,'B2');
    xlswrite(xlsfilename,z(inew),xlspage,'C2');
    xlswrite(xlsfilename,density',xlspage,'D2');
    
   clear whos density*
    clear N1 N2 N3  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Fig4: Local Density Type 3
      
    [d3_all d3_1 d3_2]=NNanalysis(x(S==3),y(S==3),z(S==3));
    
    
    inew=find(S==3);
    for i=1:length(x(S==3));       
        index=find(d3_all(i,:)<LocalRadius);%local density in a LocalRadius µm area
        density(i)=mean(d3_1(index));
    end
    
    figure
    set(gcf,'Position',[100 100 1000 800]);
    for i=1:length(x(S==3))
        ii=1+floor(63-(density(i)-ColorRange(4,1))/(ColorRange(4,2)-ColorRange(4,1))*63);
        if ii<1
            ii=1;
        elseif ii>64
            ii=64;
        end
        colorm=[CC(ii,1) CC(ii,2) CC(ii,3)];
        plot(x(inew(i)),y(inew(i)),'o','LineWidth',1,'MarkerSize',1+ii/64*4,'MarkerFaceColor',colorm,'MarkerEdgeColor',colorm);hold on;
    end
    axis equal;
    axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
    title(['Local density III (local average distance to nearest neighbor) range=[',num2str(ColorRange(4,1)),'µm ',num2str(ColorRange(4,2)),'µm] mean=',num2str(mean(d3_1),2),'µm'],'FontSize',10)
    set(gca,'LineWidth',1.5,'FontSize',10,'FontWeight','bold');    
    saveas(gcf,[path,filename(1:end-4),'Case',num2str(k),'_Analysis02Fig4'], 'tiffn');
       
    %Write to excel file
    xlspage='Fig4';
    xlswrite(xlsfilename, {'x','y','z','Local Density Type III'}, xlspage, 'A1');
    xlswrite(xlsfilename,x(inew),xlspage,'A2');
    xlswrite(xlsfilename,y(inew),xlspage,'B2');
    xlswrite(xlsfilename,z(inew),xlspage,'C2');
    xlswrite(xlsfilename,density',xlspage,'D2');
    
    
    
    
    
    
    
    
    
    
    
    
      

