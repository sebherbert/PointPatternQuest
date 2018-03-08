
function pointDynPatternMain()
% Calls and handles the output of the pointDynPatternAnalysis function.

clear


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function parameters
PARAMS = {};
PARAMS.version = 'version0p0p1';

% Datafile parameters
PARAMS.dataFile.input = uipickfiles('Prompt',...
    'Select the correct file to analyse (example: liveDataDmCurated.mat)');

% Population parameters
PARAMS.pops.pop1 = 'mother';
PARAMS.pops.pop2 = 'mother';
PARAMS.pops.popPermut = 'type2';

% Analysis parameters
PARAMS.anaMap.doMap = 1;

PARAMS.anaSearch.doMinSearch = 0;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load file table of all cells % expecting only 1 file for the moment
% File is a n*8 table of the (PositionX, PositionY,PositionZ, Unit, Time, ID, cellType, oldID)

for exp = 1:numel(PARAMS.dataFile.input)
    [PARAMS.dataFile.path{exp},PARAMS.dataFile.name{exp},~] = fileparts(PARAMS.dataFile.input{exp});
    PARAMS.dataFile.path{exp} = [PARAMS.dataFile.path{exp}, filesep];
    fullDataLive = load(PARAMS.dataFile.input{exp});
    
    % Rename data file to always be fullDataLive
    fooFN = fieldnames(fullDataLive);
    fullDataLive = fullDataLive.(fooFN{1});
    
    %% Launch analysis of the file
    pointDynPatternAnalysis()
    
    
    
    
end

end