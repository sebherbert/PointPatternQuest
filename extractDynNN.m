

function outNNdistances = extractDynNN(PARAMS,popSource,popTarget,popPermut,RSpairs)
% For one R and S pair calculate the NN distances between source population
% and the target population. For the moment, the tpSource == tptarget is
% excluded to avoid picking into the same population
%{
INPUT
PARAMs - See pointDynPatternMain
popSource - table containing 3D position of the source population
popTarget - table containing 3D position of the target population
popPermut - IF you want to simulate a model, then this variable has to be a
   table containing 3D position of the permutation population. IF this
   variable is empty then no simulation will be made

OUTPUT
NNdistances - list of NxM distances with:
   - N = 1 if experimental populations or the number of simulations to 
   effectuate (nPerms).
   - M = the number of lines in popSource

%}

switch nargin
    case 3
        if PARAMS.verbose > 1
            fprintf('Calculating NN for experimental data\n');
        end
        doSimulation = false;
    case 5 % additional info for RSpairs and popPermut
        if PARAMS.verbose > 1 
            fprintf('Calculating NN for simulated data\n');
        end
        doSimulation = true;
    otherwise
        msg = fprintf('Error: Unrecognised number of parameters in extractDynNN function.\n');
        error(msg)
end

% for every source temporal step => tpSource
for tpSource = PARAMS.movie.minTp:PARAMS.movie.maxTp-1
    % For the moment, discard the situation where tpSource = tpTarget
    shortSource = [popSource.PositionX(popSource.Time==tpSource, :),...
        popSource.PositionY(popSource.Time==tpSource, :),...
        popSource.PositionZ(popSource.Time==tpSource, :)];

    % for every target temporal step => tpTarget
    for tpTarget = tpSource+1:PARAMS.movie.maxTp
        % Change for tpTarget = tpSource:PARAMS.movie.maxTp if same tp authorized
        
        % Create a field name for data output
        tpField = sprintf('dt%d_to_dt%d',tpSource,tpTarget);
        
        % If simulation
        if doSimulation
            % Put both effect parameters into an independant variable
            effect.Range = PARAMS.anaMap.RSpairs.Rs(RSpairs);
            effect.Strength = PARAMS.anaMap.RSpairs.Ss(RSpairs);
            simulatedObjects = simulateDynDispersion(PARAMS,popSource(popSource.Time==tpSource, :),...
                popTarget(popTarget.Time==tpTarget, :),...
                popPermut(popPermut.Time==tpTarget, :),effect); % List of positions
            %             function [dnSimu, Grand] = simulateSpatialDisp(effect, NNExp,
            %             pops, rowPermut, nTarget, PARAMS) => static version
          
            allNNsimuDist = zeros(PARAMS.anaGlobal.numPermut,size(shortSource,1));
            for permN = 1:PARAMS.anaGlobal.numPermut
                shortTarget = [popPermut.PositionX(sum(popPermut.ID==simulatedObjects(permN,:),2,'native'),:),...
                    popPermut.PositionY(sum(popPermut.ID==simulatedObjects(permN,:),2,'native'),:),...
                    popPermut.PositionZ(sum(popPermut.ID==simulatedObjects(permN,:),2,'native'),:)];
                
                [~,fooNNdist] = knnsearch(shortTarget, shortSource, 'K', 1);
                allNNsimuDist(permN,:) = fooNNdist; % NN distances for this set of positions
            end
            outNNdistances.(tpField) = allNNsimuDist';
            
        else % if using the experimental population
            shortTarget = [popTarget.PositionX(popTarget.Time==tpTarget),...
                popTarget.PositionY(popTarget.Time==tpTarget),...
                popTarget.PositionZ(popTarget.Time==tpTarget)];
            
            if tpSource == tpTarget % 
                [~,fooNNdist] = knnsearch(shortTarget, shortSource, 'K', 2);
                fooNNdist = fooNNdist(:,2); % first column is identical position...
            else
                [~,fooNNdist] = knnsearch(shortTarget, shortSource, 'K', 1);
            end
            outNNdistances.(tpField) = fooNNdist;
        end % end measurement
        clear foo*
    end % end tpTarget
end % end tpSource

% 

end

%%%%%%%%%%% ADDITIONNAL FUNCTIONS %%%%%%%%%%%
% function 











