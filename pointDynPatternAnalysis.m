
function [NNdistances, NNdispersion, NNana] = pointDynPatternAnalysis(PARAMS,popSource,popTarget,popPermut)

% Test a temporal correlation in the appearance of a cell type compared to
% the appearance of an other or similar cell type.

% Associate experimental NNdistances
NNdistances.exp = extractDynNN(PARAMS,popSource,popTarget); % Doesn't depend on the RS pair
NNdispersion.exp = mergeDynNN(PARAMS, NNdistances.exp);

% Keep track of which pops are used as targets
permutPopNames = cell2mat(unique(popPermut.cellType)');


NNana = {};

%% Do correlation map
if PARAMS.anaMap.doMap
    
    % Reformat R and S settings into pairscdf
    PARAMS.anaMap.RSpairs = reformatRSpairs(PARAMS.anaMap.RminmaxnSteps, PARAMS.anaMap.SminmaxnSteps,...
        PARAMS.anaMap.Rlog, PARAMS.anaMap.Slog);   
    
    for RSpair = 1:numel(PARAMS.anaMap.RSpairs.Rs) % each R and S

        % Prepare output structure for the NNdistances with appropriate fieldnames
        pairString = sprintf('R%0.1f_S%0.3f',PARAMS.anaMap.RSpairs.Rs(RSpair), PARAMS.anaMap.RSpairs.Ss(RSpair));
        pairString = regexprep(pairString,'\.','p');
        
        if PARAMS.verbose > 0
            fprintf('Calculating NN for %s\n', pairString);
        end

        % Call Simulator for each pairs asked to simulate the expectedcell dispersion with
        % specified RSpair
        NNdistances.simus.(pairString) = extractDynNN(PARAMS, popSource, popTarget, popPermut, RSpair);
        
        % Provide a formatted metric for comparing the different simulations and
        % experimental results.
        NNdispersion.simu.(pairString) = mergeDynNN(PARAMS, NNdistances.simus.(pairString));
        
        % further analysis?
        % Use an adapted metric for the map assessment: RMSE still?
        NNana.perRSpair.(pairString) = anaDynNN(PARAMS,NNdispersion.exp,NNdispersion.simu.(pairString));
        
        % Reproduce the original histogram from WS
        if ((PARAMS.display.reproHistoWS.do) && (PARAMS.anaMap.RSpairs.Ss(RSpair) == 1))
            reproHistoWS(NNdispersion.exp, NNdispersion.simu.(pairString), PARAMS.display.deltaTOIs, permutPopNames);
            PARAMS.display.reproHistoWS.do = 0;
        end
        
    end
    
    toc

    % Display the populations
    %     NNdisplays = {}; % Initialize structure
    DTfields = fieldnames(NNdispersion.exp); % Max number of Δt in the movie (frame-1)

    PARAMS.anaMap.uniqR = unique(PARAMS.anaMap.RSpairs.Rs);
    PARAMS.anaMap.uniqS = unique(PARAMS.anaMap.RSpairs.Ss);

    for DTfield = PARAMS.display.deltaTOIs(1):PARAMS.display.deltaTOIs(end) % interesting Δts
        % Reparse the different pairs of R and S for the display (to keep things separated...)        
                
        if PARAMS.anaGlobal.doRMSE 
            % Reshape and add further analyses per deltaT
            NNana.perDeltaT.(DTfields{DTfield}).RMSE = reshapeMap(PARAMS, ...
                NNana, 'RMSE', PARAMS.anaMap.uniqS, PARAMS.anaMap.uniqR, (DTfields{DTfield}));
            statTest = 'RMSE';
        end
    end
    
    if PARAMS.display.NNmap % Launch display of maps
        displayDynNNmap(PARAMS, PARAMS.anaMap.uniqR, PARAMS.anaMap.uniqS,...
            DTfields, statTest, permutPopNames, NNana.perDeltaT);
    end
    if PARAMS.display.NNisoMap.do % Launch display of iso maps
        displayDynNNisomap(PARAMS, PARAMS.anaMap.uniqR, PARAMS.anaMap.uniqS,...
            DTfields, statTest, permutPopNames, NNana.perDeltaT);
    end
    
        
    save(sprintf('%s/NNdistances',PARAMS.dataFile.currentPath),'NNdistances');
    save(sprintf('%s/NNdispersion',PARAMS.dataFile.currentPath),'NNdispersion');
    save(sprintf('%s/NNanalysis',PARAMS.dataFile.currentPath),'NNana');    
end


%% Do search function
if PARAMS.anaSearch.doMinSearch
    % Call Simulator search tool
    fprintf('Section under construction...');
end


end



function outNNana = reshapeMap(PARAMS, NNana, statTest, Ss, Rs, DTfieldName)
% Adds a structure to the analysis organized per deltaT instead of RS pair
% the statistical test represents an individual field

outNNana = {};

RSfieldNames = fieldnames(NNana.perRSpair);
for RSpair = 1:numel(PARAMS.anaMap.RSpairs.Rs) % each R and S

    outNNana.map(RSpair) = NNana.perRSpair.(RSfieldNames{RSpair}).(statTest).(DTfieldName).(statTest);
    
end % endfor RSpair

outNNana.map = reshape(outNNana.map, numel(Ss), numel(Rs));
[outNNana.minSloc, outNNana.minRloc] = find(outNNana.map == min(outNNana.map(:)));
outNNana.minSval = Ss(outNNana.minSloc);
outNNana.minRval = Rs(outNNana.minRloc);
outNNana.minVal = outNNana.map(outNNana.minSloc, outNNana.minRloc);

end



function reproHistoWS(inNNdispersionExp, inNNdispersionSimu, deltaTOIs, permutPopNames)
% This function aims a the reproduction of the already provided simu
% histogram vs exp average when no effect are set in the model.
% display only the deltaTOIs deltaT. (array of double)

figure

numSubP = numSubplots(numel(deltaTOIs));

for deltaN = 1:numel(deltaTOIs)
    subplot(numSubP(1),numSubP(2),deltaN)
    
    deltaTOI = sprintf('deltaT%d',deltaTOIs(deltaN)); % delta Time Of Interest

    avExp = mean(inNNdispersionExp.(deltaTOI).allDistances);
    hold on;
    histogram(inNNdispersionSimu.(deltaTOI).allDistances);
    line([avExp, avExp], ylim, 'LineWidth', 2, 'Color', 'r');
    
    title(sprintf('NN for deltat=%dtp', deltaTOIs(deltaN)))
    xlabel('Distance (µm)'); ylabel('# sim N');
    
    legend({'simulated data' sprintf('exp av = %0.1fµm',avExp)})
end

saveas(gcf,sprintf('reproHistWS_%s).fig',permutPopNames));

end


function RSpairs = reformatRSpairs(RminmaxnSteps, SminmaxnSteps, Rlog, Slog)
% reformat R and S settings into pairs for linear list handling

if Rlog % use logarithm spacing (base 10)
    allR = logspace(log10(RminmaxnSteps(1)),log10(RminmaxnSteps(2)),RminmaxnSteps(3));
else % use linear spacing
    allR = linspace(RminmaxnSteps(1),RminmaxnSteps(2),RminmaxnSteps(3));
end

if Slog % use logarithm spacing (base 10)
    allS = logspace(log10(SminmaxnSteps(1)),log10(SminmaxnSteps(2)),SminmaxnSteps(3));
else % use linear spacing
    allS = linspace(SminmaxnSteps(1),SminmaxnSteps(2),SminmaxnSteps(3));
end

RSpairs = table;
for Ri = 1:numel(allR)
    for Si = 1:numel(allS)
        RSpairs = vertcat(RSpairs, table(allR(Ri),allS(Si)));
    end
end

RSpairs.Properties.VariableNames = {'Rs', 'Ss'};

end



