% Source is DataAnalysis.m from Willy Supatto modified Sept 27, 2017
% Script copied and modified by Sebastien Herbert


function pointPatternMain()

clear
close all

PARAMS = {};
PARAMS.dispDistrib_1 = 0;
PARAMS.dispDensityMap_2 = 0;
PARAMS.numPermut = 200;

PARAMS.effect = 'None'; 
% Type of effect of a cell on its nearest neighbours can only be 'None',
% 'Repulsion', 'Attraction'
PARAMS.range = 1; % Distance of effect of a cell on its neighbours
PARAMS.strength = 2; % Strength of the effect of a cell on its neighbours

% Display Parameters
PARAMS.maxSizeCDF = 200; % maximum number of points on the cdf
PARAMS.binSize = 0:0.1:14; % bin size for the ecdf => if force binning of ecdf

% File import
fileToOpen = uipickfiles('Prompt','Please, select the correct file to analyse (example: sox2_C_subdiv_L_corrected_nodb_noDl.ims)','num',1);
% for development purposes only
% fileToOpen = {'/media/sherbert/Data/Projects/OG_projects/Project6_ND_distribPattern/static_/preAnalysis_and_ims/sox2_C_subdiv_L_corrected_nodb_noDl_xyzCASE1.mat'};
[path,name,ext] = fileparts(fileToOpen{1});
path = [path, filesep];
filename = [name,ext];

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

% find which part of the brain is studied
% regexp => single digit following the '_xyzCASE' string
k = str2double(regexp(name, '(?<=_xyzCASE)[0-9]','match'));
switch k
    case 1 
        brainPart = '_allSample';
    case 2
        brainPart = '_da';
    case 3
        brainPart = '_dl';
    case 4
        brainPart = '_dm';
    otherwise
        warning('Unexpected part of he brain. Stopping the analysis');
        return
end
fprintf('Analysing the part of the brain: %s (case %d)\n',brainPart,k);

dataCombined = {};
dataCombined.name = name;
dataCombined.brainPart = brainPart;


if exist([path,name,ext]) % useless in the end... 
    disp(name);
    load([path,name,ext]);
    
    % S = cell type
    % dN_M = distance in population N between cell of interest and Mth
    % neighbour
    % 
    
    if PARAMS.dispDistrib_1
        %% ANALYSIS 01: Distribution Display
        %Parameters Analysis 01
        ColorRange(1)=7;%Fig 3&4: min distance in �m (red)
        ColorRange(2)=12;%Fig 3&4: max distance in �m (blue
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Main Function
        Analysis01DistributionDisplay(path,filename,k,x,y,z,S,mean(d123_1),mean(d12_1),mean(d13_1),ColorRange)
        clear ColorRange
    end
    
    if PARAMS.dispDensityMap_2
        %% ANALYSIS 02: Density Maps
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
        %Note: LocalRadius=50 t        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %         oo small in this case...
        
        %Local averaging radius (50 �m)
        LocalRadius=50;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Main Function
        Analysis02DensityMaps(path,filename,k,x,y,z,S,d12_all,d123_all,d12_1,d123_1,LocalRadius,ColorRange)
        clear ColorRange
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Point Pattern Analysis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Main Functions
    
    %% These functions are not ready yet but the rest needed to be merged
    %     NNExp = table((1:length(S))',S',d123_1',[x,y,z]);
    %     NNExp.Properties.VariableNames = {'cellID','cellType','nearestNeighbour','pos3D'};
    %
    %     % Point pattern analysis Type 2 effect on Type 2 in Type 1+2 (first neighbor)
    %     pops.popSource = 2;
    %     pops.popTarget = 2;
    %     pops.popPermut = [1 2];
    %     tAnalysis = sprintf('t%dvst%d',pops.popSource,pops.popTarget);
    %     fullPath = [path, name, '_', tAnalysis];
    %     fprintf('\nRunning analysis %s\n',tAnalysis);
    %     dataCombined.(tAnalysis) = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS);
    %
    %     % Point pattern analysis Type 3 effect on type 3 in Type 1+2+3 (first neighbor)
    %     pops.popSource = 3;
    %     pops.popTarget = 3;
    %     pops.popPermut = [1 2 3];
    %     tAnalysis = sprintf('t%dvst%d',pops.popSource,pops.popTarget);
    %     fullPath = [path, name, '_', tAnalysis];
    %     fprintf('\nRunning analysis %s\n',tAnalysis);
    %     dataCombined.(tAnalysis) = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS);
    %
    %     % Point pattern analysis Type 3 effect on type 2 in Type 1+2 (first neighbor)
    %     pops.popSource = 3;
    %     pops.popTarget = 2;
    %     pops.popPermut = [1 2];
    %     tAnalysis = sprintf('t%dvst%d',pops.popSource,pops.popTarget);
    %     fullPath = [path, name, '_', tAnalysis];
    %     fprintf('\nRunning analysis %s\n',tAnalysis);
    %     dataCombined.(tAnalysis) = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS);
    
    dataCombined.t2vst2 = analysis03PPAFirstNeighbor(path,name,k,x,y,z,S,d12_all,d12_1,PARAMS);
    %Point pattern analysis Type 2 in Type 1+2 (first neighbor)
    
    %         Analysis04PointPatternAnalysisSecondNeighbor(path,filename,k,x,y,z,S,d12_all,d12_1)
    %Point pattern analysis Type 2 in Type 1+2 (second neighbor)
    
    dataCombined.t3vst3 = analysis05PPAFirstNeighbor(path,name,k,x,y,z,S,d123_all,d123_1,PARAMS);
    %Point pattern analysis Type 3 in Type 1+2+3 (first neighbor)
    
    %     Analysis06PointPatternAnalysisSecondNeighbor(path,filename,k,x,y,z,S,d123_all,d123_1)
    %Point pattern analysis Type 3 in Type 1+2+3 (second neighbor)
    
    dataCombined.t3vst2 = analysis07PPAFirstNeighbor(path,name,k,x,y,z,S,d123_all,d123_1,PARAMS);
    %Point pattern analysis Type 2 distance to 3 in Type 1+2 (first neighbor)
    
end

save([path,name,brainPart],'dataCombined');

end
