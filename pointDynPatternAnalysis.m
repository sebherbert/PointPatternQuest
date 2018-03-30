
function [NNdistances, NNdispersion]= pointDynPatternAnalysis(PARAMS,popSource,popTarget,popPermut)

% Test a temporal correlation in the appearance of a cell type compared to
% the appearance of an other or similar cell type.

% Associate experimental NNdistances
NNdistances.exp = extractDynNN(PARAMS,popSource,popTarget); % Doesn't depend on the RS pair
NNdispersion.exp = analyseDynNN(PARAMS, NNdistances.exp);

if PARAMS.anaMap.doMap == 1
    
    % Reformat R and S settings into pairs
    PARAMS.anaMap.RSpairs = reformatRSpairs(PARAMS.anaMap.RminmaxnSteps, PARAMS.anaMap.SminmaxnSteps,...
        PARAMS.anaMap.Rlog, PARAMS.anaMap.Slog);   
    
    for RSpair = 1:numel(PARAMS.anaMap.RSpairs.Rs) % each R and S

        % Prepare output structure for the NNdistances with appropriate fieldnames
        pairString = sprintf('R%0.1f_S%0.2f',PARAMS.anaMap.RSpairs.Rs(RSpair), PARAMS.anaMap.RSpairs.Ss(RSpair));
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
        
        % Reproduce the original histogram from WS
        
        %         reproHistoWS()
        
        % Use an adapted metric for the map assessment: RMSE still?
        
        % foo = reshape(inNNdispersionSimu.deltaT4.allDistances,[],1);
        % figure
        % cdfplot(foo)
        % hold on
        % cdfplot(inNNdispersionExp.deltaT4.allDistances)

        
        
        
    end
    
    toc
    
    save(sprintf('%s/NNdistances',PARAMS.dataFile.path{1}),'NNdistances');
    save(sprintf('%s/NNdispersion',PARAMS.dataFile.path{1}),'NNdispersion');
end



if PARAMS.anaSearch.doMinSearch == 1
    % Call Simulator search tool
    fprintf('Section under construction...');
end


end


%
% function reproHistoWS(inNNdispersionExp, inNNdispersionSimu, deltaTOIs)
% % This function aims a the reproduction of the already provided simu
% % histogram vs exp average when no effect are set in the model.
% % display only the deltaTOIs deltaT. (array of double)
%
% deltaTOI = 'deltaT1'; % delta Time Of Interest
%
%
% figure
% avExp = mean(inNNdispersionExp.(deltaTOI).allDistances);
% hold on;
% foo = histogram(inNNdispersionSimu.(deltaTOI).allDistances);
% line([avExp, avExp], ylim, 'LineWidth', 2, 'Color', 'r');
%
% end

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



