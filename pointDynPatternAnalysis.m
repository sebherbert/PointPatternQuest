
function [NNdistances, NNdispersion, NNana, PARAMS] = pointDynPatternAnalysis(PARAMS,popSource,popTarget,popPermut)

% Test a temporal correlation in the appearance of a cell type compared to
% the appearance of an other or similar cell type.
% Should separate from display but no time right now...

% Keep track of which pops are used as targets
permutPopNames = cell2mat(unique(popPermut.cellType)');

%% Do correlation map
if PARAMS.anaMap.doMap
    
    if PARAMS.anaMap.loadMap % if the map is already calculated. Requires the exact same parameters!
        if exist('completeAnalysis.mat')==2 % check if there is already a completeAnalysis file in the current folder
            fooLoad = load('completeAnalysis.mat');
        else % if not, asks for a new input
            file2Load = uipickfiles('Prompt',...
                'Select the data file to display (example: completeAnalysis.mat)', 'NumFile', 1);
            fooLoad = load(file2Load{1}); % Can only be one file
        end
        fooFN = fieldnames(fooLoad);
        
        if strcmp(fooFN{1},'completeAnalysis') % Check the data structure
            NNdistances = fooLoad.(fooFN{1}).NNdistances;
            NNdispersion = fooLoad.(fooFN{1}).NNdispersion;
            fooLoad.(fooFN{1}).NNPARAMS.dummy = PARAMS.dummy; % force to update the dummy parameter to match the new call
            PARAMS = fooLoad.(fooFN{1}).NNPARAMS;
        else
            error('Expecting a "completeAnalysis" first field, check that the right file was selected. (pointDynPatternmain.m output)');
        end
        
        % Clearing temp var
        clear foo*;

    else
        % Reformat R and S settings into pairscdf
        PARAMS.anaMap.RSpairs = reformatRSpairs(PARAMS.anaMap.RminmaxnSteps, PARAMS.anaMap.SminmaxnSteps,...
            PARAMS.anaMap.Rlog, PARAMS.anaMap.Slog);
        
        % Associate experimental NNdistances
        NNdistances.exp = extractDynNN(PARAMS,popSource,popTarget); % Doesn't depend on the RS pair
        NNdispersion.exp = mergeDynNN(PARAMS, NNdistances.exp);
        
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
                       
            %             % Reproduce the original histogram from WS
            %             if ((PARAMS.display.reproHistoWS.do) && (PARAMS.anaMap.RSpairs.Ss(RSpair) == 1))
            %                 reproHistoWS(PARAMS, NNdispersion.exp, NNdispersion.simu.(pairString), ...
            %                     PARAMS.display.deltaTOIs, permutPopNames);
            %                 PARAMS.display.reproHistoWS.do = 0;
            %             end
            
        end
    end % end of loadMap test
    
    toc
    
    NNana = {};
    fieldNs = fieldnames(NNdispersion.simu);
    for RSpair = 1:numel(fieldNs) % each R and S do the analysis
        % Use an adapted metric for the map assessment: RMSE still?
        NNana.perRSpair.(fieldNs{RSpair}) = ...
            anaDynNN(PARAMS,NNdispersion.exp,NNdispersion.simu.(fieldNs{RSpair}));
    end
    
    toc
    
    % Display the populations
    DTfields = fieldnames(NNdispersion.exp); % Max number of Δt in the movie (frame-1)
    
    PARAMS.anaMap.uniqR = unique(PARAMS.anaMap.RSpairs.Rs);
    PARAMS.anaMap.uniqS = unique(PARAMS.anaMap.RSpairs.Ss);

    for DTfield = PARAMS.display.deltaTOIs(1):PARAMS.display.deltaTOIs(end) % interesting Δts
        
        % Reparse the different pairs of R and S for the display (to keep things separated...)
        NNana.perDeltaT.(DTfields{DTfield}).exp = NNexpAna(NNdispersion.exp.(DTfields{DTfield}));
        
        if PARAMS.anaGlobal.doRMSE
            % Reshape and add further analyses per deltaT
            statTest = 'RMSE';
            NNana.perDeltaT.(DTfields{DTfield}).RMSE = reshapeMap(PARAMS, ...
                NNana, statTest, PARAMS.anaMap.uniqS, PARAMS.anaMap.uniqR, (DTfields{DTfield}));
        end
    end
    
    % Reproduce the original histogram from WS
    if PARAMS.display.reproHistoWS.do
        reproHistoWS(PARAMS, NNdispersion, NNana, DTfields, permutPopNames);
    end    
    
    if PARAMS.display.NNmap % Launch display of maps
        iso = 0; % don't show isocontours
        displayDynNNmap(PARAMS, PARAMS.anaMap.uniqR, PARAMS.anaMap.uniqS,...
            DTfields, statTest, permutPopNames, NNana.perDeltaT, iso);
    end
    if PARAMS.display.NNisoMap.do % Launch display of iso maps
        iso = 1; % show isocontours
        displayDynNNmap(PARAMS, PARAMS.anaMap.uniqR, PARAMS.anaMap.uniqS,...
            DTfields, statTest, permutPopNames, NNana.perDeltaT, iso);
    end
     
    if ~PARAMS.dummy
    %     if PARAMS.anaMap.saveOutputMat
        save(sprintf('%s/NNdistances',PARAMS.dataFile.currentPath),'NNdistances');
        save(sprintf('%s/NNdispersion',PARAMS.dataFile.currentPath),'NNdispersion');
        save(sprintf('%s/NNanalysis',PARAMS.dataFile.currentPath),'NNana');
    end
end


%% Do search function
if PARAMS.anaSearch.doMinSearch
    % Call Simulator search tool
    fprintf('Section under construction...');
end


end


function outNNana = NNexpAna(expDists)
% Extract some simple analyses from the experimetal values

outNNana.Ncells = numel(expDists.allDistances);
% average based measurements
outNNana.average = mean(expDists.allDistances);
outNNana.std = std(expDists.allDistances);
outNNana.sem = outNNana.std / sqrt(outNNana.Ncells);
% median based measurements
outNNana.median = median(expDists.allDistances);
outNNana.iqr = iqr(expDists.allDistances);

end

function outNNana = reshapeMap(PARAMS, inNNana, statTest, Ss, Rs, DTfieldName)
% Adds a structure to the analysis organized per deltaT instead of RS pair
% the statistical test represents an individual field

RSfieldNames = fieldnames(inNNana.perRSpair);
for RSpair = 1:numel(PARAMS.anaMap.RSpairs.Rs) % each R and S

    outNNana.map(RSpair) = inNNana.perRSpair.(RSfieldNames{RSpair}).(statTest).(DTfieldName).(statTest);
    
end % endfor RSpair

outNNana.map = reshape(outNNana.map, numel(Ss), numel(Rs));
[outNNana.minSloc, outNNana.minRloc] = find(outNNana.map == min(outNNana.map(:)));
outNNana.minSval = Ss(outNNana.minSloc);
outNNana.minRval = Rs(outNNana.minRloc);
outNNana.minVal = outNNana.map(outNNana.minSloc, outNNana.minRloc);
end



function reproHistoWS(PARAMS, NNdispersion, NNana, DTfields, permutPopNames)
% This function aims a the reproduction of the already provided simu
% histogram vs exp average when no effect are set in the model.
% display only the deltaTOIs deltaT. (array of double)
% It has the possibility to merge all the strength model = 1 (No effect)

mergedModels = {};

allSimus = fieldnames(NNdispersion.simu);
for simu = 1:numel(allSimus) % for each simulation
    modelStrength = str2double(regexprep(regexprep(allSimus{simu},'\w*_S',''),'p','.'));
    if modelStrength == 1 % If the strength of the model is 1 (no effect)
        for deltaN = 1:numel(DTfields) % fill each deltaN 
            if isfield(mergedModels, DTfields{deltaN}) % create field if it doens't exist or concatenate it
                mergedModels.(DTfields{deltaN}) = horzcat(mergedModels.(DTfields{deltaN}), ...
                    NNdispersion.simu.(allSimus{simu}).(DTfields{deltaN}).allDistances);
            else
                mergedModels.(DTfields{deltaN}) = NNdispersion.simu.(allSimus{simu}).(DTfields{deltaN}).allDistances;
            end
        end
    end
end

%% Display average model against average exp
figure
hold on;
numSubP = numSubplots(numel(PARAMS.display.deltaTOIs)); % prepare subplot values
for deltaN = 1:numel(PARAMS.display.deltaTOIs)
    subplot(numSubP(1),numSubP(2),deltaN)   
    simuData = mean(mergedModels.(DTfields{deltaN}),1); % Calculate average distance in the model
    avVsHisto(simuData, DTfields, NNana, deltaN); % Display the image for deltaT#
end

if ~PARAMS.dummy
    saveas(gcf,sprintf('averageDistFreq_%s.fig',permutPopNames));
end

%% Display all distances of model against average exp
figure
hold on;
numSubP = numSubplots(numel(PARAMS.display.deltaTOIs)); % prepare subplot values
for deltaN = 1:numel(PARAMS.display.deltaTOIs)
    subplot(numSubP(1),numSubP(2),deltaN)   
    simuData = reshape(mergedModels.(DTfields{deltaN}),[],1); % Calculate all distances in the model
    avVsHisto(simuData, DTfields, NNana, deltaN); % Display the image for deltaT#
end

if ~PARAMS.dummy
    saveas(gcf,sprintf('allDistFreq_%s.fig',permutPopNames));
end

% 
% figure
% hold on
% ecdf(simuData(:,1))
% ecdf(simuData(:,2))
% ecdf(simuData(:,3))


end

function avVsHisto(simuData, DTfields, NNana, deltaN)

% h{1} = histogram(simuData,'BinWidth',5,'Normalization','pdf','EdgeColor','None');
h{1} = histogram(simuData,'BinMethod','integers','Normalization','pdf','EdgeColor','none');
h{1}.NumBins = 100;

% display experimental results
avExp = NNana.perDeltaT.(DTfields{deltaN}).exp.average;
semExp = NNana.perDeltaT.(DTfields{deltaN}).exp.sem;
x = [avExp-semExp, avExp+semExp, avExp+semExp, avExp-semExp];
y = [0 0 max(h{1}.Values) max(h{1}.Values)];
h{2} = line([avExp, avExp], [0 max(h{1}.Values)],  ylim, 'Color', 'r','LineStyle','-');
h{3} = line([avExp-semExp, avExp+semExp ; avExp-semExp, avExp+semExp], [0 0 ;max(h{1}.Values) max(h{1}.Values)],  ylim, 'Color', 'r','LineStyle','--');
patch(x, y, 'r', 'EdgeColor','None','FaceAlpha',0.1)

uistack(h{1},'top');
title(sprintf('NN for %s (N=%d)', DTfields{deltaN}, NNana.perDeltaT.(DTfields{deltaN}).exp.Ncells));
xlabel('Distance (µm)'); ylabel('<Dist. Freq>');

legend([h{1} h{2} h{3}(1,:)], {'simulated data' sprintf('exp av = %0.1fµm',avExp) sprintf('exp sem +/- %0.1fµm',semExp)})

end


function RSpairs = reformatRSpairs(RminmaxnSteps, SminmaxnSteps, Rlog, Slog)
% reformat R and S settings into pairs for linear list handling

if Rlog % use logarithm spacing (base 10)
    allR = logspace(log10(RminmaxnSteps(1)),log10(RminmaxnSteps(2)),RminmaxnSteps(3));
else % use linear spacing
    allR = linspace(RminmaxnSteps(1),RminmaxnSteps(2),RminmaxnSteps(3));
end

if Slog % use logarithm spacing (base 10)
    allS = logspace(log10(SminmaxnSteps(1)),log10(SminmaxnSteps(2)),SminmaxnSteps(3))';
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



