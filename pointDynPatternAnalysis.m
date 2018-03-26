
function analysisOut = pointDynPatternAnalysis(PARAMS,popSource,popTarget,popPermut)

% Test a temporal correlation in the appearance of a cell type compared to
% the appearance of an other or similar cell type.

% Associate experimental NNdistances
NNdistances.exp = movieAnalysisNN(PARAMS,popSource,popTarget); % Doesn't depend on the RS pair

if PARAMS.anaMap.doMap == 1
    
    % Reformat R and S settings into pairs
    PARAMS.anaMap.RSpairs = reformatRSpairs(PARAMS.anaMap.RminmaxnSteps, PARAMS.anaMap.SminmaxnSteps,...
        PARAMS.anaMap.Rlog, PARAMS.anaMap.Slog);   
    
    for RSpair = 1:numel(PARAMS.anaMap.RSpairs.Rs) % each R and S

        % Prepare output structure for the NNdistances with appropriate fieldnames
        pairString = sprintf('R%0.1f_S%0.2f',PARAMS.anaMap.RSpairs.Rs(RSpair), PARAMS.anaMap.RSpairs.Ss(RSpair));
        pairString = regexprep(pairString,'\.','p');

        % Call Simulator for each pairs asked to simulate the movie analysis with
        % specified RSpair
        NNdistances.simus.(pairString) = movieAnalysisNN(PARAMS, popSource, popTarget, popPermut, RSpair);
        
        % Create CDFs out of NNdistances ; format into proper shape
        
        
        % create RMSE map
        %         outMapAna =
    end
end



if PARAMS.anaSearch.doMinSearch == 1
    % Call Simulator search tool
    fprintf('Section under construction...');
end


end



function calculationRMSEMap()

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



