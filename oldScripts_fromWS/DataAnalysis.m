%Modified Sept 27, 2017 

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TO CHECK
path = 'D:\Spatial Analysis Nico Laure\temp nico\';
filename = 'sox2_C_subdiv_L_corrected_nodb_noDl.ims';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%REMINDER
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

%%
for k=1
     %  for    k=4
    if exist([path,filename(1:end-4),'_xyzCASE',num2str(k),'.mat'])
        [path,filename(1:end-4),'_xyzCASE',num2str(k),'.mat']
        load([path,filename(1:end-4),'_xyzCASE',num2str(k),'.mat']);

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANALYSIS 01: Distribution Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        %Parameters Analysis 01
        ColorRange(1)=7;%Fig 3&4: min distance in �m (red)
        ColorRange(2)=12;%Fig 3&4: max distance in �m (blue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Main Function
        Analysis01DistributionDisplay(path,filename,k,x,y,z,S,mean(d123_1),mean(d12_1),mean(d13_1),ColorRange)
        clear ColorRange
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANALYSIS 02: Density Maps
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Parameters Analysis 02
        %Fig1: local density I+II+III
        ColorRange(1,1)=5;%Fig 1: min distance in �m (red)
        ColorRange(1,2)=9;%Fig 1: max distance in �m (blue)
        
        %Fig2: local density I+II      
        ColorRange(2,1)=5;%Fig 2: min distance in �m (red)
        ColorRange(2,2)=9;%Fig 2: max distance in �m (blue)
                      
        %Fig3: local density II
        ColorRange(3,1)=10;%Fig 3: min distance in �m (red)
        ColorRange(3,2)=50;%Fig 3: max distance in �m (blue)
        
        %Fig4: local density III
        ColorRange(4,1)=7;%Fig 4: min distance in �m (red)
        ColorRange(4,2)=20;%Fig 4: max distance in �m (blue)
        
        %Fig5: local ratio II/I+II+III
        ColorRange(5,1)=0;%Fig 5: min % (blue)
        ColorRange(5,2)=12;%Fig 5: max % (red)
        
        %Fig6: local ratio III/I+II+III
        ColorRange(6,1)=0;%Fig 6: min % (blue)
        ColorRange(6,2)=20;%Fig 6: max % (red)
        
        %Fig7: local ratio II+III/I+II+III 
        ColorRange(7,1)=0;%Fig 7: min % (blue)
        ColorRange(7,2)=30;%Fig 7: max % (red)
        
        %Fig8: local ratio II/I+II        
        ColorRange(8,1)=0;%Fig 8: min % (blue)
        ColorRange(8,2)=12;%Fig 8: max % (red)
        
        %Fig9: local ratio II/I+II        
        ColorRange(9,1)=0;%Fig 9: min % (blue)
        ColorRange(9,2)=100;%Fig 9: max % (red)        
        %Note: LocalRadius=50 too small in this case...
       
        %Local averaging radius (50 �m)
        LocalRadius=50;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Main Function
        Analysis02DensityMaps(path,filename,k,x,y,z,S,d12_all,d123_all,d12_1,d123_1,LocalRadius,ColorRange)
        clear ColorRange
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%ANALYSIS 03 - 06: Point Pattern Analysis
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Main Functions

        Analysis03PointPatternAnalysisFirstNeighbor(path,filename,k,x,y,z,S,d12_all,d12_1)
        %Point pattern analysis Type 2 in Type 1+2 (first neighbor)

        Analysis04PointPatternAnalysisSecondNeighbor(path,filename,k,x,y,z,S,d12_all,d12_1)
        %Point pattern analysis Type 2 in Type 1+2 (second neighbor)
        
        Analysis05PointPatternAnalysisFirstNeighbor(path,filename,k,x,y,z,S,d123_all,d123_1)
        %Point pattern analysis Type 3 in Type 1+2+3 (first neighbor)
        
        Analysis06PointPatternAnalysisSecondNeighbor(path,filename,k,x,y,z,S,d123_all,d123_1)
        %Point pattern analysis Type 3 in Type 1+2+3 (second neighbor)

        Analysis07PointPatternAnalysisFirstNeighbor(path,filename,k,x,y,z,S,d123_all,d123_1)
        %Point pattern analysis Type 2 distance to 3 in Type 1+2 (first neighbor)
      
    end
    
    
end

close all
