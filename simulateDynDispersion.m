
function outSimObjects = simulateDynDispersion(PARAMS,popSource,popTarget,popPermut,effect)
% Simulates the expected dispersion for a population among an other
% population based on the parameters Range and Strength of the effect.
% WARNING! Is NOT meant to sustain a large amount of targets, will 
% otherwise consume too much memory! 
%{
INPUT:
PARAMS - Structure of relevant parameters
popSource - Table of cells at the origin of an effect
popTarget - Table of experimental cells affected by the Source cells
popPermut - Table of cells among which new popTarget will be drawn

OUTPUT:
NNdistances - list of NxM objects with:
   - N = the number of simulations to effectuate (nPerms).
   - M = the number of lines in popTarget

%}

nTarget = height(popTarget);

outSimObjects = zeros(PARAMS.anaGlobal.numPermut , height(popTarget));

% Adapt the probability map once and for all (as long as diff pop)
probMap = table();
probMap.proba = ones(height(popPermut),1);
probMap.ID = cellstr(num2str(popPermut.ID));

probMap = adaptProbMap(popSource, popPermut, effect, probMap);

popPermut.proba = probMap.proba;

for permutation = 1:PARAMS.anaGlobal.nPermut % Will become a parfor at some point!
    
    if mod(perm,500) == 0
        fprintf('Running permutation %d\n',perm);
    end
    
    simulatedPopTargetIDs = simulateDispersionIDs(effect, popSource, popPermut, nTarget, PARAMS); % List of simulated population IDs
    outSimObjects(permutation,:) = simulatedPopTargetIDs;
end

end

function simulatedPopTargetIDs = simulateDispersionIDs(effect, popSource, popPermut, nTarget, PARAMS)
% Simulates random draws of the popTarget population among the popPermut
% population, considering the effect of the popSource population on the
% next draw.

IDsSimu = zeros(1,nTarget);

% Initiate a probability map => prob = 1 for the permutable population and
% 0 for the rest of the population
popTargetSimu{perm} = table;
%
%     if PARAMS.samePop => As long as the time point can not be the same,
%     it can not be the same population.
%         % Requires adaptation of the popSource at every step + an adaptation at
%         % every step of the draw
%         for bioCell = 1:nTarget
%             tempCell = datasample(NNExp,1,'Weights',probMap);
%             % Take cell out of the permutable population
%             % be careful at not putting it back at more than 0 by an attractive or repulsive effect !
%             probMap(tempCell.cellID) = 0;
%
%             % Adapt the probability map of the other cells based on the draw.
%             probMap = adaptProbMap(tempCell.cellID, probMap, cell2CellDist, effect);
%
%
%             % Add the drawn cell to the complete simulated draw
%             popTargetSimu{perm} = [popTargetSimu{perm};tempCell];
%         end
%         % In this case target and source populations are the same
%         popSourceSimu{perm} = popTargetSimu{perm};
%
%     else

%         
%         % Then make a single draw without replacing the cells
%         popTargetSimu{perm} = datasample(NNExp,nTarget,'Replace',false,'Weights',probMap);
%         popSourceSimu{perm} = NNExp(NNExp.cellType == pops.popSource, :);
%     end
%     
%     % Find the Nearest neighbours in the population
%     dnSimu(:,perm) = findNN(popSourceSimu{perm}.pos3D, popTargetSimu{perm}.pos3D, PARAMS.samePop, pops);
%     
%     %     % Calculate the cdf histogram
%     %     for i=1:size(PARAMS.binSize,2)
%     %         Grand(perm,i)= sum(dnSimu(:,perm)<PARAMS.binSize(i)) / numel(dnSimu(:,perm));
%     %     end
% end

end

function newProbMap = adaptProbMap(popSource, popPermut, effect, oriProbMap)
% Adapt the draw probability map based on the desired effect, the drawned
% cell(s)
% popSource = table of the source population affecting the probability map
% popPermut = table of the affected permutation population
% effect of the cell type (.Range / . Strength)

tempProbMap = oriProbMap;

% Find objects at a distance < effect.Range
listID = rangesearch([popSource.PositionX, popSource.PositionY, popSource.PositionZ],...
    [popPermut.PositionX, popPermut.PositionY, popPermut.PositionZ], effect.Range);
% listID: popPermut lines each stating which of the popSource are in range

% Change their probability
for bioCell = 1:height(popPermut)
    for sourceInRange = 1:numel(listID{bioCell})
        tempProbMap.proba(bioCell) = tempProbMap.proba(bioCell)*effect.Strength;
    end
end

% if later choose to change the multiplication factor by an additive term,
% multiply tempProMap by a vector containing the 0 positions of the
% oriProbMap (1 otherwise) to ensure conservation of the forbidden
% positions => only for same data analysis
newProbMap = tempProbMap;

end