%Update sep 29, 2017

clear all
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO CHECK
path = 'D:\Spatial Analysis Nico Laure\temp nico\';
filename = 'sox2_C_subdiv_L_corrected_nodb_noDl.ims';


! "C:\Program Files\Bitplane\Imaris x64 9.0.2\Imaris.exe" id101 &
% javaaddpath('C:\Program Files\Bitplane\Imaris x64 9.0.2\XT\rtmatlab\ImarisLib.jar');
javaaddpath('C:\Program Files\Bitplane\Imaris x64 9.0.2\XT\rtmatlab\ImarisLib.jar');


%Export to .mat file?
export=1;%0=NO 1=YES

%Save .tif figures?
saveFig=1;%0=NO 1=YES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



pause(6);

ImarisId = 101;

vImarisLib = ImarisLib;
vImarisApplication = vImarisLib.GetApplication(ImarisId);
vImarisApplication.FileOpen([path,filename],'');

%%

vScene = vImarisApplication.GetSurpassScene;
vFactory = vImarisApplication.GetFactory;
vNumberOfChildren = vScene.GetNumberOfChildren;

for vChildIndex = 1:vNumberOfChildren
    
    vObject = vScene.GetChild(vChildIndex-1);
    
    % We check if the child is a Reference Frame
    if strcmp(vObject.GetName,'BodyAxes')
        
        referenceFrame = vImarisApplication.GetFactory.ToReferenceFrames(vScene.GetChild(vChildIndex-1));
        quatStruct = referenceFrame.GetFramesQuaternion(0);
        quat = quatStruct.mQuaternions;
        q0 = quat(4);
        q1 = -quat(1);
        q2 = -quat(2);
        q3 = -quat(3);
        
        newX = [(1-2*q2*q2-2*q3*q3) 2*(q1*q2-q0*q3) 2*(q1*q3+q0*q2)];
        newY = [2*(q1*q2+q0*q3) (1-2*q1*q1-2*q3*q3) 2*(q3*q2-q0*q1)];
        newZ = [2*(q3*q1-q0*q2) 2*(q3*q2+q0*q1) (1-2*q1*q1-2*q2*q2)];
        
        pa = newX/norm(newX);
        rl = newY/norm(newY);
        vd = newZ/norm(newZ);
        
    end
    
    % We check if the child is a DataContainer
    if vFactory.IsDataContainer(vScene.GetChild(vChildIndex-1))
        Type = vFactory.ToDataContainer(vObject);
        
        if strcmp(vObject.GetName, 'type 1');
            type1_sox2p = vFactory.ToSpots(Type.GetChild(0)).GetPositionsXYZ;
            type1_sox2m = vFactory.ToSpots(Type.GetChild(1)).GetPositionsXYZ;
            type1_sox2p_da = vFactory.ToSpots(Type.GetChild(2)).GetPositionsXYZ;
            type1_sox2m_da = vFactory.ToSpots(Type.GetChild(3)).GetPositionsXYZ;
            type1_sox2p_dl = vFactory.ToSpots(Type.GetChild(4)).GetPositionsXYZ;
            type1_sox2m_dl = vFactory.ToSpots(Type.GetChild(5)).GetPositionsXYZ;
            type1_sox2p_dm = vFactory.ToSpots(Type.GetChild(6)).GetPositionsXYZ;
            type1_sox2m_dm = vFactory.ToSpots(Type.GetChild(7)).GetPositionsXYZ;
        elseif strcmp(vObject.GetName, 'type 2');
            type2_sox2p = vFactory.ToSpots(Type.GetChild(0)).GetPositionsXYZ;
            type2_sox2m = vFactory.ToSpots(Type.GetChild(1)).GetPositionsXYZ;
            type2_sox2p_da = vFactory.ToSpots(Type.GetChild(2)).GetPositionsXYZ;
            type2_sox2m_da = vFactory.ToSpots(Type.GetChild(3)).GetPositionsXYZ;
            type2_sox2p_dl = vFactory.ToSpots(Type.GetChild(4)).GetPositionsXYZ;
            type2_sox2m_dl = vFactory.ToSpots(Type.GetChild(5)).GetPositionsXYZ;
            type2_sox2p_dm = vFactory.ToSpots(Type.GetChild(6)).GetPositionsXYZ;
            type2_sox2m_dm = vFactory.ToSpots(Type.GetChild(7)).GetPositionsXYZ;
        elseif strcmp(vObject.GetName, 'type 3');
            type3_sox2p = vFactory.ToSpots(Type.GetChild(0)).GetPositionsXYZ;
            type3_sox2m = vFactory.ToSpots(Type.GetChild(1)).GetPositionsXYZ;
            type3_sox2p_da = vFactory.ToSpots(Type.GetChild(2)).GetPositionsXYZ;
            type3_sox2m_da = vFactory.ToSpots(Type.GetChild(3)).GetPositionsXYZ;
            type3_sox2p_dl = vFactory.ToSpots(Type.GetChild(4)).GetPositionsXYZ;
            type3_sox2m_dl = vFactory.ToSpots(Type.GetChild(5)).GetPositionsXYZ;
            type3_sox2p_dm = vFactory.ToSpots(Type.GetChild(6)).GetPositionsXYZ;
            type3_sox2m_dm = vFactory.ToSpots(Type.GetChild(7)).GetPositionsXYZ;
        end
        
        
    end
    
end


% save([path,filename(1:end-4),'_xyzDATA.mat'],'-regexp', 'type1','-regexp', 'type2','-regexp', 'type3','pa','vd','rl');
%%



a=who('-regexp', 'type');
for i=1:length(a)
     eval(['position=',a{i},';']);
     if size(position,1)==1&&sum(position==[0 0 0])==3
         eval([a{i},'=[];']);
     end
end


data(1).type(1).sox2m=type1_sox2m;
data(1).type(1).sox2p=type1_sox2p;
data(1).type(2).sox2m=type2_sox2m;
data(1).type(2).sox2p=type2_sox2p;
data(1).type(3).sox2m=type3_sox2m;
data(1).type(3).sox2p=type3_sox2p;
data(2).type(1).sox2m=type1_sox2m_da;
data(2).type(1).sox2p=type1_sox2p_da;
data(2).type(2).sox2m=type2_sox2m_da;
data(2).type(2).sox2p=type2_sox2p_da;
data(2).type(3).sox2m=type3_sox2m_da;
data(2).type(3).sox2p=type3_sox2p_da;
data(3).type(1).sox2m=type1_sox2m_dl;
data(3).type(1).sox2p=type1_sox2p_dl;
data(3).type(2).sox2m=type2_sox2m_dl;
data(3).type(2).sox2p=type2_sox2p_dl;
data(3).type(3).sox2m=type3_sox2m_dl;
data(3).type(3).sox2p=type3_sox2p_dl;
data(4).type(1).sox2m=type1_sox2m_dm;
data(4).type(1).sox2p=type1_sox2p_dm;
data(4).type(2).sox2m=type2_sox2m_dm;
data(4).type(2).sox2p=type2_sox2p_dm;
data(4).type(3).sox2m=type3_sox2m_dm;
data(4).type(3).sox2p=type3_sox2p_dm;

%changement de referentiel
pa=pa';rl=rl';vd=vd';
% pa=[0 1 0]';rl=[1 0 0]';vd=[0 0 1]';
P=[rl pa vd];
pa=P\pa;
rl=P\rl;
vd=P\vd;

%%
for k=1:4
    
    %all:k=1 da:k=2 dl:k=3 dm:k=4
    sample=[data(k).type(1).sox2m;data(k).type(1).sox2p;data(k).type(2).sox2m;data(k).type(2).sox2p;data(k).type(3).sox2m;data(k).type(3).sox2p];
    
    
    if ~isempty(sample)
        
        x=sample(:,1);y=sample(:,2);z=sample(:,3);
        
        %%%%%
        N1m=size(data(k).type(1).sox2m,1);
        N1p=size(data(k).type(1).sox2p,1);
        N2m=size(data(k).type(2).sox2m,1);
        N2p=size(data(k).type(2).sox2p,1);
        N3m=size(data(k).type(3).sox2m,1);
        N3p=size(data(k).type(3).sox2p,1);
        
        Type=[];
        Type(1:N1m)=10;
        Type(N1m+1:N1m+N1p)=11;
        Type(N1m+N1p+1:N1m+N1p+N2m)=20;
        Type(N1m+N1p+N2m+1:N1m+N1p+N2m+N2p)=21;
        Type(N1m+N1p+N2m+N2p+1:N1m+N1p+N2m+N2p+N3m)=30;
        Type(N1m+N1p+N2m+N2p+N3m+1:N1m+N1p+N2m+N2p+N3m+N3p)=31;
        
        S=Type;
        S(Type==10|Type==11)=1;
        S(Type==20|Type==21)=2;
        S(Type==30|Type==31)=3;
        
        %Type
        %10=type 1 sox2-
        %11=type 1 sox2+
        %20=type 2 sox2-
        %21=type 2 sox2+
        %30=type 3 sox2-
        %31=type 2 sox2+
        
        %S
        %1=type 1 (or Type 10 + Type 11)
        %2=type 2 (or Type 20 + Type 31)
        %3=type 3 (or Type 30 + Type 31)
        
        %CASE k
        % k=1 all
        % k=2 da
        % k=3 dl
        % k=4 dm
        
        
        
        
        %Check duplicated points
        A=sample;
        Nspots=size(A,1);
        Duplicata=x.*0;
        id=x'.*0;
        
        clear d
        
        for i=1:Nspots
            id=id.*0;id(i)=1;
            if Duplicata(i)==0
                for j=1:Nspots;
                    d(j)=norm([x(i)-x(j) y(i)-y(j) z(i)-z(j)]);
                end
                
                if length(find(d==0))>2%in case we have triplets or more
                    find(d==0)
                end
                Duplicata(find(d==0&id==0))=1;
            end
        end
        
        if sum(Duplicata)>0
            warning(['case #',num2str(k),' ',num2str(sum(Duplicata)),' duplicated cells!']);
        end
        
        
        
        
        
        
        
        
        
        
        
        
        figure(k)
        % plot(y,x,'.k');hold on
        
        for i=1:size(x,1)
            B=P\[x(i,1);y(i,1);z(i,1)];
            x(i,1)=B(1);
            y(i,1)=-B(2);
            z(i,1)=B(3);
        end
        
        xMin(k)=min(x);
        yMin(k)=min(y);
        zMin(k)=min(z);
        
        x=x-xMin(1);
        y=y-yMin(1);
        z=z-zMin(1);
        plot(x(S==1),y(S==1),'o','LineWidth',1,'MarkerSize',1,'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0]);hold on;
        plot(x(S==2),y(S==2),'o','LineWidth',2,'MarkerSize',2,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0]);hold on;
        plot(x(S==3),y(S==3),'o','LineWidth',2,'MarkerSize',2,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 1 0]);hold on;
        
        axis equal;
        axis([min(x)-20 max(x)+20 min(y)-20 max(y)+20])
        
        title(['case #',num2str(k),' N_{I+II+III}=',num2str(size(x,1)),' N_{II}=',num2str(length(S(S==2))),' (',num2str(100*length(S(S==2))/size(x,1),2),'% red)',' N_{III}=',num2str(length(S(S==3))),' (',num2str(100*length(S(S==3))/size(x,1),2),'% green)'])
        
        if saveFig==1
            saveas(gcf,[path,filename(1:end-4),'_DataCheckFig',num2str(k)], 'tiffn');
        end
        
        if export==1
            save([path,filename(1:end-4),'_xyzCASE',num2str(k),'.mat'],'x','y','z','Type','S');
        end
        
        
    end
end

if export==1
    for k=1:4
        if exist([path,filename(1:end-4),'_xyzCASE',num2str(k),'.mat'])
            load([path,filename(1:end-4),'_xyzCASE',num2str(k),'.mat']);
            
            %Type1+Type2+Type3
            [d123_all,d123_1,d123_2]=NNanalysis(x,y,z);
            %Type1+Type2
            [d12_all,d12_1,d12_2]=NNanalysis(x(S==1|S==2),y(S==1|S==2),z(S==1|S==2));
            %Type1+Type3
            [d13_all,d13_1,d13_2]=NNanalysis(x(S==1|S==3),y(S==1|S==3),z(S==1|S==3));
            %         %Type2 only
            %         [d2_all d2_1 d2_2]=NNanalysis(x(S==2),y(S==2),z(S==2));
            
            save([path,filename(1:end-4),'_xyzCASE',num2str(k),'.mat'], '-regexp', 'd','-append');
        end
    end
end

close all 