
function pointDynPatternMain()
% Calls and handles the output of the pointDynPatternAnalysis function.

clear


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function parameters
PARAMS = {};
PARAMS.version = 'version0p0p1';

PARAMS.verbose = 1;

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
PARAMS.anaGlobal.numPermut = 200; % Number of random permutations
% RMSE map
PARAMS.anaMap.doMap = 1; % => RMSE map
PARAMS.anaMap.RminmaxnSteps = [5 30 11]; % R min value, max value, number of steps
PARAMS.anaMap.SminmaxnSteps = [1/16 16 17]; % S min value, max value, number of steps
PARAMS.anaMap.Rlog = false; % If you want to have a logarithmic scale in R
PARAMS.anaMap.Slog = true; % If you want to have a logarithmic scale in S
% 
PARAMS.anaSearch.doMinSearch = 0; % => min search protocol

% Movie parameters
PARAMS.movie.maxTp = NaN; % To be filled up after loading the movie
PARAMS.movie.minTp = NaN;
PARAMS.movie.dt = NaN; % To be filled up after  



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

    % Call analysis launcher for each independant file
    completeAnalysis = PDPACaller(PARAMS, fullDataLive);
    
    % Save PDPA
    save(completeAnalysis);
end
     
    
end


function completeAnalysis = PDPACaller(PARAMS, fullDataLive)
% handles the parsing of the different populations in the different rounds
% of analysis before sending them to the main analysis function

for anaSeq = 1:height(PARAMS.pops) % => for each analysis sequence

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
    for permuts = 1:size(PARAMS.pops.popPermut{anaSeq},1)
        popPermut = vertcat(popPermut, fullDataLive(...
            strcmp(fullDataLive.cellType,PARAMS.pops.popPermut{permuts}),:));
    end
    
    % Launch the pattern analysis
    analysisOut = pointDynPatternAnalysis(PARAMS,popSource,popTarget,popPermut);
    
    % Displays
    
    % Merge to previous analyses
    
end
    


end










