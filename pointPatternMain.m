% Source is DataAnalysis.m from Willy Supatto modified Sept 27, 2017
% Script copied and modified by Sebastien Herbert


function pointPatternMain()

clear

PARAMS = {};
PARAMS.version = 'version0p1p3';
PARAMS.dispDistrib_1 = 0;
PARAMS.dispDensityMap_2 = 0;
PARAMS.numPermut = 500;

% Optimization parameters
PARAMS.optimizePar = 1; % Do an automated search for the best parameters
PARAMS.minFitRange = 5; % Minimum Range for the fitted model
PARAMS.minFitStrength = 0; % Minimum Strength for the fitted model
% Original values for the optimization function
PARAMS.optiR0 = 10; % 10µm distance
PARAMS.optiS0 = 1; % No dispersion effect
PARAMS.fitMaxIter = 50;

PARAMS.displayIndivCDF = 0; % => To display individual cdf vs model figures

% Type of effect of a cell on its nearest neighbours can only be 'None',
% 'Repulsion', 'Attraction'
% Distance of effect of a cell on its neighbours
% PARAMS.effectRange = 10;
PARAMS.effectMultiRange = [10];

PARAMS.effectRangeU = 'µm';

% Strength of the effect of a cell on its neighbours
% x = inputdlg('Enter space-separated numbers:',...
%     'Sample', [1 50]);
% PARAMS.effectMultiStrength = str2num(x{:}); % Can be multiple values
PARAMS.effectMultiStrength = [1]; % Can be multiple values
% if effect is > 1 => it increases the probability of the neighbours to be
% selected => Attraction
% if effect is < 1 => it decreases the probability of the neighbours to be
% selected => Repulsion
% if effect is = 1 => it doesn't affect the probability of the neighbours to be
% selected => No effect



% Display Parameters
PARAMS.maxSizeCDF = 200; % maximum number of points on the cdf
PARAMS.binSize = 0:0.1:100; % bin size for the ecdf => if force binning of ecdf
PARAMS.axis = [0 100 0 1];

% File import
fileToOpen = uipickfiles('Prompt','Please, select the correct file to analyse (example: sox2_C_subdiv_L_corrected_nodb_noDl_xyzCASE1.mat)');
% for development purposes only
% fileToOpen = {'/media/sherbert/Data/Projects/OG_projects/Project6_ND_distribPattern/static_/preAnalysis_and_ims/sox2_C_subdiv_L_corrected_nodb_noDl_xyzCASE1.mat'};

tic

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

dataCombinedModels = {};
dataCombinedModels.name = PARAMS.name;
dataCombinedModels.brainPart = PARAMS.brainPart;


disp(PARAMS.name);
load([PARAMS.path,PARAMS.name,ext]);

% Begin the PPA analysis for each model
if PARAMS.optimizePar % aka if you want the fitted version
    % Change the model type into a model name for output
    PARAMS.model  = 'fittedModel';
    dataCombinedModels.fittedModel = mainPPA(S, d123_1, x, y, z, PARAMS);

    dataCombined = dataCombinedModels.(PARAMS.model);
    
    save([PARAMS.path,PARAMS.name,PARAMS.brainPart,'_',PARAMS.model],'dataCombined');

else
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
            dataCombinedModels.(PARAMS.model) = mainPPA(S, d123_1, x, y, z, PARAMS);
            
            dataCombinedModels.(PARAMS.model).name = PARAMS.name;
            dataCombinedModels.(PARAMS.model).brainPart = PARAMS.brainPart;
            dataCombinedModels.(PARAMS.model).effectType = PARAMS.effectType;
            dataCombinedModels.(PARAMS.model).effectStrength = PARAMS.effectStrength;
            dataCombinedModels.(PARAMS.model).effectRange = PARAMS.effectRange;
            dataCombinedModels.(PARAMS.model).effectRangeU = PARAMS.effectRangeU;
            
            dataCombined = dataCombinedModels.(PARAMS.model);
            
            save([PARAMS.path,PARAMS.name,PARAMS.brainPart,'_',PARAMS.model],'dataCombined');            
        end
        
    end
    
    if numel(PARAMS.effectMultiStrength)*numel(PARAMS.effectMultiRange) > 1
        displayModelSerie(dataCombinedModels, PARAMS);
        save([PARAMS.path,PARAMS.name,PARAMS.brainPart,'_allModels'],'dataCombinedModels');
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

% Point pattern analysis Type 2 effect on Type 2 in Type 1+2 (first neighbor)
pops.popSource = 2;
pops.popTarget = 2;
pops.popPermut = [1 2];
tAnalysis = sprintf('t%dvst%d',pops.popSource,pops.popTarget);
fullPath = [PARAMS.path, PARAMS.name, '_', tAnalysis];
fprintf('Running analysis %s\n',tAnalysis);
dataCombined.(tAnalysis) = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS);

% Point pattern analysis Type 3 effect on type 3 in Type 1+2+3 (first neighbor)
pops.popSource = 3;
pops.popTarget = 3;
pops.popPermut = [1 2 3];
tAnalysis = sprintf('t%dvst%d',pops.popSource,pops.popTarget);
fullPath = [PARAMS.path, PARAMS.name, '_', tAnalysis];
fprintf('Running analysis %s\n',tAnalysis);
dataCombined.(tAnalysis) = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS);

% Point pattern analysis Type 3 effect on type 2 in Type 1+2 (first neighbor)
pops.popSource = 3;
pops.popTarget = 2;
pops.popPermut = [1 2];
tAnalysis = sprintf('t%dvst%d',pops.popSource,pops.popTarget);
fullPath = [PARAMS.path, PARAMS.name, '_', tAnalysis];
fprintf('Running analysis %s\n',tAnalysis);
dataCombined.(tAnalysis) = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS);

end
