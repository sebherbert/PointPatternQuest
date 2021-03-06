% Source is DataAnalysis.m from Willy Supatto modified Sept 27, 2017
% Script copied and modified by Sebastien Herbert
% Model: 
% if Strength is > 1 => it increases the probability of the neighbours to be
% selected => Attraction
% if Strength is < 1 => it decreases the probability of the neighbours to be
% selected => Repulsion
% if Strength is = 1 => it doesn't affect the probability of the neighbours to be
% selected => No effect

function pointPatternMain()

clear


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS/ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PARAMS = {};
PARAMS.version = 'version0p1p11';

%% Which cell type distribution effect to test
PARAMS.dot2vst2 = 1;
PARAMS.dot3vst3 = 0;
PARAMS.dot3vst2 = 1;

%% Global Display parameters
PARAMS.dispDistrib_1 = 0;
PARAMS.dispDensityMap_2 = 0;
PARAMS.displayIndivCDF = 0; % => To display individual cdf vs model figures
PARAMS.maxSizeCDF = 200; % maximum number of points on the cdf
PARAMS.binSize = 0:0.1:100; % bin size for the ecdf => if force binning of ecdf in µm
PARAMS.axis = [0 40 0 1];
PARAMS.dispModelsOverlay = 0; % When different Range or Strength are tested 

%% Global saving parameters
PARAMS.saveIndivModel = 0; % When different Range or Strength are tested
PARAMS.suffix = '_allTests'; % add a suffix to the filename of the save

%% RMSE calculation parameters
PARAMS.numPermut = 1000; % Nbr of permutation for model estimation % Default is 1000
PARAMS.useRMSMaxDist = 1; % if you want to use only a part of the RMS (as a function of the NN)
PARAMS.maxDistFactor = 3; % times the cellDiameter for the RMS limit.
if PARAMS.useRMSMaxDist % Automatically updates suffix name
    PARAMS.suffix = sprintf('%s_max%dcellDia', PARAMS.suffix, PARAMS.maxDistFactor);
end
PARAMS.cellDiameter = nan; % Will be set in Analysis function

%% Optimization model parameters (fitting method)
PARAMS.doOptimFit = 1; % Do an automated search for the best parameters
PARAMS.minFitRange = 5; % Minimum Range for the fitted model
PARAMS.minFitStrength = 0; % Minimum Strength for the fitted model
% Original values for the optimization function
PARAMS.doVarFitInit = 0; % => optiR0 and optiS0 will be changing the folder name
PARAMS.doDisplayLiveFit = 0; % If you want to see the fit evolve live
PARAMS.useRangeCDF50 = 1; % Boolean to exchange the cdf 
PARAMS.optiR0 = 10; % µm Range => WARNING WILL BE OVERWRITTEN IF doVarFitInit
PARAMS.optiS0 = 1; % Dispersion Strength  => WARNING WILL BE OVERWRITTEN IF doVarFitInit
PARAMS.fitMaxIter = 200; % Nbr of iterations for the min search % Default is 200

%% Hard coded model parameters (map method)
% Type of effect of a cell on its nearest neighbours can only be 'None',
% 'Repulsion', 'Attraction'

PARAMS.doMapRMSE = 0; % Do a map of specific models
% Distance of effect of a cell on its neighbours
PARAMS.effectMultiRange = 5.1:0.2:30; % Can be multiple values => Default
% PARAMS.effectMultiRange = 5.1:2:20; %  => for dev values
PARAMS.effectRangeU = 'µm';
% Strength of the effect of a cell on its neighbours
PARAMS.effectMultiStrength = logspace(log(1/16)/log(10),log(16)/log(10),17); % Can be multiple values
% PARAMS.effectMultiStrength = 1; % Can be multiple values

%% Merge cells that are too close together
PARAMS.cellMerge.do = 0; % Do or do not
PARAMS.cellMerge.cellType = {2}; % Which cell type?
PARAMS.cellMerge.distThresh = 5; % Distance below which cell pairs should be merged

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% /PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% File import
fileToOpen = uipickfiles('Prompt','Please, select the correct file to analyse (example: sox2_C_subdiv_L_corrected_nodb_noDl_xyzCASE1.mat)');
% for development purposes only
% fileToOpen = {'/media/sherbert/Data/Projects/OG_projects/Project6_ND_distribPattern/static_/preAnalysis_and_ims/sox2_C_subdiv_L_corrected_nodb_noDl_xyzCASE1.mat'};

tic

%% Apply function to individual files
for fileOfInterest = 1:length(fileToOpen)
    mainBodyFunctions(fileToOpen{fileOfInterest},PARAMS);
end

toc

end

function mainBodyFunctions(dataFile,PARAMS)

close all

[PARAMS.path,PARAMS.name,ext] = fileparts(dataFile);
PARAMS.path = [PARAMS.path, filesep];
filename = [PARAMS.name,ext];

cd(PARAMS.path)

if (PARAMS.doVarFitInit && PARAMS.doOptimFit) % adapts fit initialization based on folder name
    folderIndexes = strfind(PARAMS.path,filesep) ;
    lastFolderName = PARAMS.path(folderIndexes(end-1)+1:folderIndexes(end)-1);
    if contains(lastFolderName,{'R','S'})
        PARAMS.optiR0 = str2num(cell2mat(regexp(lastFolderName, '(?<=R)[0-9]*','match')));
        if strcmp(cell2mat(regexp(lastFolderName, '(?<=S)[0-9]*','match')),'0')
            PARAMS.optiS0 = str2num(cell2mat(regexprep(regexp(lastFolderName, '(?<=S)[0-9]*0p[0-9]*','match'),'p','.')));
        else
            PARAMS.optiS0 = str2num(cell2mat(regexp(lastFolderName, '(?<=S)[0-9]*','match')));
        end
        fprintf('Settings new fit init values to R0 = %d and S0 = %d\n', PARAMS.optiR0, PARAMS.optiS0);
    else 
        fprintf('folder name %s is inadapted! Stopping the analysis\n',lastFolderName);
        return 
    end
end

% find which part of the brain is studied
% regexp => single digit following the '_xyzCASE' string
k = str2double(regexp(PARAMS.name, '(?<=_xyzCASE)[0-9]','match'));
switch k
    case 1
        PARAMS.brainPart = '_allSample';
    case 2
        PARAMS.brainPart = '_da';
    case 3
        PARAMS.brainPart = '_dl';
    case 4
        PARAMS.brainPart = '_dm';
    otherwise
        warning('Unexpected part of the brain. Stopping the analysis');
        return
end
fprintf('Analysing the part of the brain: %s (case %d)\n',PARAMS.brainPart,k);


disp(PARAMS.name);
load([PARAMS.path,PARAMS.name,ext]);

% Begin the PPA analysis for each model
% if PARAMS.optimizePar % aka if you want the fitted version
if PARAMS.doOptimFit % aka if you want the fitted version
    % Change the model type into a model name for output
    PARAMS.model  = 'fittedModel';

    dataCombined_fitModel.data = mainPPA(S, d123_1, x, y, z, PARAMS);
    dataCombined_fitModel.PARAMS = PARAMS;
       
    save([PARAMS.path,PARAMS.name,PARAMS.brainPart,PARAMS.suffix,'_',PARAMS.model],'dataCombined_fitModel');
end

if PARAMS.doMapRMSE
    
    dataCombined_RMSMap = {};
    dataCombined_RMSMap.name = PARAMS.name;
    dataCombined_RMSMap.brainPart = PARAMS.brainPart;
    PARAMS.model  = 'RMSmap';
    dataCombined_RMSMap.PARAMS = PARAMS;
    
    for modelR = 1:numel(PARAMS.effectMultiRange)
        % For every range tested
        PARAMS.effectRange = PARAMS.effectMultiRange(modelR);
        
        for modelS = 1:numel(PARAMS.effectMultiStrength)
            % For every strength tested
            PARAMS.effectStrength = PARAMS.effectMultiStrength(modelS);
            
            if PARAMS.effectStrength > 1
                PARAMS.effectType = 'Attraction';
            elseif PARAMS.effectStrength < 1
                PARAMS.effectType = 'Repulsion';
            else
                PARAMS.effectType = 'None';
            end
            
            % Change the model type into a model name for output
            tempRange = regexprep(num2str(PARAMS.effectRange),'\.','p');
            tempStrength = regexprep(num2str(PARAMS.effectStrength),'\.','p');
            PARAMS.model = sprintf('Model_T%s_R%s_S%s',PARAMS.effectType(1),...
                tempRange,tempStrength);
            fprintf('\nProcessing model %s\n',PARAMS.model);
            dataCombined_RMSMap.(PARAMS.model) = mainPPA(S, d123_1, x, y, z, PARAMS);
            
            dataCombined_RMSMap.(PARAMS.model).name = PARAMS.name;
            dataCombined_RMSMap.(PARAMS.model).brainPart = PARAMS.brainPart;
            dataCombined_RMSMap.(PARAMS.model).effectType = PARAMS.effectType;
            dataCombined_RMSMap.(PARAMS.model).effectStrength = PARAMS.effectStrength;
            dataCombined_RMSMap.(PARAMS.model).effectRange = PARAMS.effectRange;
            dataCombined_RMSMap.(PARAMS.model).effectRangeU = PARAMS.effectRangeU;
            
            dataCombined = dataCombined_RMSMap.(PARAMS.model);
            
            if PARAMS.saveIndivModel
                save([PARAMS.path,PARAMS.name,PARAMS.brainPart,'_',PARAMS.model],'dataCombined');
            end
        end
        
    end
    
    if numel(PARAMS.effectMultiStrength)*numel(PARAMS.effectMultiRange) > 1
        if PARAMS.dispModelsOverlay
            displayModelSerie(dataCombined_RMSMap, PARAMS);
        end
        
        save([PARAMS.path,PARAMS.name,PARAMS.brainPart,PARAMS.suffix, '_RMSMap'],'dataCombined_RMSMap');
        % (=> still gonna give 3analType*3brainPart / sample)
        % => Adapt Combine and compare ? Only the display part ?
        % => Just copy and modify it ? keep the KS tests ?
    end
end

end

function dataCombined = mainPPA(S, d123_1, x, y, z, PARAMS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Point Pattern Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Main Functions

NNExp = table((1:length(S))',S',d123_1',[x,y,z]);
NNExp.Properties.VariableNames = {'cellID','cellType','nearestNeighbour','pos3D'};

% Merge cells that are too close together
if PARAMS.cellMerge.do
    NNExp = mergeCloseCells(PARAMS, NNExp);
end

% Point pattern analysis Type 2 effect on Type 2 in Type 1+2 (first neighbor)
if PARAMS.dot2vst2
    pops.popSource = 2;
    pops.popTarget = 2;
    pops.popPermut = [1 2];
    tAnalysis = sprintf('t%dvst%d',pops.popSource,pops.popTarget);
    fullPath = [PARAMS.path, PARAMS.name, '_', tAnalysis];
    fprintf('Running analysis %s\n',tAnalysis);
    dataCombined.(tAnalysis) = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS);
end

% Point pattern analysis Type 3 effect on type 3 in Type 1+2+3 (first neighbor)
if PARAMS.dot3vst3
    pops.popSource = 3;
    pops.popTarget = 3;
    pops.popPermut = [1 2 3];
    tAnalysis = sprintf('t%dvst%d',pops.popSource,pops.popTarget);
    fullPath = [PARAMS.path, PARAMS.name, '_', tAnalysis];
    fprintf('Running analysis %s\n',tAnalysis);
    dataCombined.(tAnalysis) = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS);
end

% Point pattern analysis Type 3 effect on type 2 in Type 1+2 (first neighbor)
if PARAMS.dot3vst2
    pops.popSource = 3;
    pops.popTarget = 2;
    pops.popPermut = [1 2];
    tAnalysis = sprintf('t%dvst%d',pops.popSource,pops.popTarget);
    fullPath = [PARAMS.path, PARAMS.name, '_', tAnalysis];
    fprintf('Running analysis %s\n',tAnalysis);
    dataCombined.(tAnalysis) = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS);
end

end


