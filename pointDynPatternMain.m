
function pointDynPatternMain()
% Calls and handles the output of the pointDynPatternAnalysis function.

clear


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function parameters
PARAMS = {};
PARAMS.version = 'version1p0p0';

PARAMS.verbose = 1;

PARAMS.dummy = 1; % do not save anything

% Datafile parameters
PARAMS.dataFile.input = uipickfiles('Prompt',...
    'Select the correct file to analyse (example: liveDataDmCurated.mat)');

% Population parameters => Add an additionnal line for each additionnal
% analysis
PARAMS.pops = table;
PARAMS.pops = vertcat(PARAMS.pops,cell2table({{'mother'}, {'mother'}, {'type2'}}));
PARAMS.pops = vertcat(PARAMS.pops,cell2table({{'mother'}, {'mother'}, {'type2';'type1'}}));
PARAMS.pops.Properties.VariableNames = {'popSource' 'popTarget' 'popPermut'};

% Analysis parameters
PARAMS.anaGlobal.numPermut = 10000; % Number of random permutations
PARAMS.anaGlobal.doRMSE = 1; % Test for RMSE between experimental and simulation curves
% Correlation map
PARAMS.anaMap.doMap = 1; % Process the correlation map
PARAMS.anaMap.loadMap = 1; % If the simulation has already been done and just needs display
PARAMS.anaMap.saveOutputMat = 0; % Save the output mat files
PARAMS.anaMap.RminmaxnSteps = [5 50 19]; % 19 R min value, max value, number of steps
PARAMS.anaMap.SminmaxnSteps = [1/64 64 13]; % 25 S min value, max value, number of steps
PARAMS.anaMap.Rlog = false; % If you want to have a logarithmic scale in R
PARAMS.anaMap.Slog = true; % If you want to have a logarithmic scale in S
% Optimization function
PARAMS.anaSearch.doMinSearch = 0; % => min search protocol

% Movie parameters
PARAMS.movie.maxTp = NaN; % To be filled up after loading the movie
PARAMS.movie.minTp = NaN; % To be filled up after loading the movie (if for some reason tp1 is ~= 1...)
PARAMS.movie.dt = NaN; % To be filled up after  

% Display parameters - map analysis
PARAMS.display.deltaTOIs = [1 2 3 4 5 6];
PARAMS.display.reproHistoWS.do = 1;
PARAMS.display.reproHistoWS.useAllPermut = 1; % Will merge all the S=1 simulations
PARAMS.display.NNmap = 1;
PARAMS.display.NNisoMap.do = 1;
PARAMS.display.NNisoMap.interFactor = 5; % If you want to interpolate



% Saving parameters
% PARAMS.saving.merge = 1; % Will overwrite if not => to improve



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

%% Load file table of all cells % expecting only 1 file for the moment
% File is a n*8 table of the (PositionX, PositionY,PositionZ, Unit, Time, ID, cellType, oldID)

for exper = 1:numel(PARAMS.dataFile.input)
    [PARAMS.dataFile.path{exper},PARAMS.dataFile.name{exper},~] = fileparts(PARAMS.dataFile.input{exper});
    PARAMS.dataFile.path{exper} = [PARAMS.dataFile.path{exper}, filesep];
    fullDataLive = load(PARAMS.dataFile.input{exper});
    
    % Rename data file to always be fullDataLive
    fooFN = fieldnames(fullDataLive);
    fullDataLive = fullDataLive.(fooFN{1});
    
    % Set PARAMS to fit the movie of interest
    PARAMS.movie.maxTp = max(fullDataLive.Time);
    PARAMS.movie.minTp = min(fullDataLive.Time); % in case of time point at 0...
    
    % Move to the folder of the movie of interest
    cd(PARAMS.dataFile.path{exper});
    PARAMS.dataFile.currentPath = PARAMS.dataFile.path{exper};

    % Call analysis launcher for each independant file
    completeAnalysis = PDPACaller(PARAMS, fullDataLive);
    
end
    
toc
    
end


function completeAnalysis = PDPACaller(PARAMS, fullDataLive)
% handles the parsing of the different populations in the different rounds
% of analysis before sending them to the main analysis function

for anaSeq = 1:height(PARAMS.pops) % => for each analysis sequence

    % Create a new dump folder
    % find out in which pops the targets will be permuted
    if size(PARAMS.pops.popPermut{anaSeq},1) > 1
        permutPopNames = strjoin(PARAMS.pops.popPermut{anaSeq});
    else
        permutPopNames = PARAMS.pops.popPermut{anaSeq}; 
    end

    % make sure you are in the right folder again...
    newFoldName = sprintf('analysisPermut%s',permutPopNames);
    mkdir(newFoldName);
    PARAMS.dataFile.currentPathRoot = PARAMS.dataFile.currentPath;
    PARAMS.dataFile.currentPath = strcat(PARAMS.dataFile.currentPath,filesep,newFoldName);
    cd(PARAMS.dataFile.currentPath);
    
    popSource = table; % fill up source table
    for sources = 1:size(PARAMS.pops.popSource{anaSeq},1)
        popSource = vertcat(popSource, fullDataLive(...
            strcmp(fullDataLive.cellType,PARAMS.pops.popSource{sources}),:));
    end

    popTarget = table; % fill up target table
    for targets = 1:size(PARAMS.pops.popTarget{anaSeq},1)
        popTarget = vertcat(popTarget, fullDataLive(...
            strcmp(fullDataLive.cellType,PARAMS.pops.popTarget{targets}),:));
    end
    
    popPermut = table; % fill up permutation table
    if size(PARAMS.pops.popPermut{anaSeq},1) == 1
        popPermut = vertcat(popPermut, fullDataLive(...
            strcmp(fullDataLive.cellType,PARAMS.pops.popPermut{anaSeq}),:));
    else
        for permuts = 1:size(PARAMS.pops.popPermut{anaSeq},1)
            popPermut = vertcat(popPermut, fullDataLive(...
                strcmp(fullDataLive.cellType,PARAMS.pops.popPermut{anaSeq}{permuts}),:));
        end
    end
    
    % Launch the pattern analysis
    [NNdistances, NNdispersion, NNana, NNPARAMS] = ...
        pointDynPatternAnalysis(PARAMS,popSource,popTarget,popPermut);
       
    % Merge to previous analyses
    completeAnalysis.NNdistances = NNdistances;
    completeAnalysis.NNdispersion = NNdispersion;
    completeAnalysis.NNana = NNana;
    completeAnalysis.NNPARAMS = NNPARAMS;
    
    % Save PDPA
    if PARAMS.anaMap.saveOutputMat
        save('completeAnalysis','completeAnalysis');
    end
    
    % Go back to root folder
    cd(PARAMS.dataFile.currentPathRoot);
    PARAMS.dataFile.currentPath = PARAMS.dataFile.currentPathRoot;    
end
    


end










