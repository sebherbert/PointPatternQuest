

function displayDynNNisomap(PARAMS, Rs, Ss, DTfields, statTest, permutPopNames, inNNana)
% all maps of a specific stat test is displayed in a dedicated subplot for
% each deltaT in an isosurface

% interpolate the 2D map to smooth it

figure

numSubP = numSubplots(numel(PARAMS.display.deltaTOIs));

% Display maps for each DeltaTs
for DTfield = PARAMS.display.deltaTOIs(1):PARAMS.display.deltaTOIs(end)

    statMap = reshape(inNNana.(DTfields{DTfield}).(statTest).map, numel(Ss), numel(Rs));

    %     % interpolate stat map using the new dimensional arrays => contour does it automatically 
    %     interpStatMap = interp2(statMap, PARAMS.display.NNisoMap.interFactor);
    %     [interRs, interSs] = interStatMap(PARAMS.anaMap.RminmaxnSteps, PARAMS.anaMap.SminmaxnSteps,...
    %         PARAMS.anaMap.Rlog, PARAMS.anaMap.Slog, ...
    %         size(interpStatMap,1), size(interpStatMap,2));

    
    subplot(numSubP(1),numSubP(2),DTfield)
    hold on
    
    contour(Rs, Ss, statMap)
    %     contour(interRs, interSs, interpStatMap);
    
    scatter3(inNNana.(DTfields{DTfield}).(statTest).minRval,...
        inNNana.(DTfields{DTfield}).(statTest).minSval,...
        inNNana.(DTfields{DTfield}).(statTest).minVal , 'filled','MarkerFaceColor',[217 83 25]/255);
    title(DTfields{DTfield});
    xlabel('Range'); ylabel('Strength');
    set(gca, 'YScale', 'log');
    legend({sprintf('%s isomap',statTest) sprintf('Global minimum = %0.3f',...
        inNNana.(DTfields{DTfield}).(statTest).minVal)},...
        'Location', 'northeast');
    colorbar

end

 saveas(gcf,sprintf('%sisomap_permutIn_%s.fig',statTest, permutPopNames));

end


function [interRs, interSs] = interStatMap(RminmaxnSteps, SminmaxnSteps, Rlog, Slog, NS, NR)
% Recalculates a new set of Rs and Ss positions for which the stat map will
% be interpolated


if Rlog % use logarithm spacing (base 10)
    interRs = logspace(log10(RminmaxnSteps(1)),log10(RminmaxnSteps(2)),NR);
else % use linear spacing
    interRs = linspace(RminmaxnSteps(1),RminmaxnSteps(2),NR);
end

if Slog % use logarithm spacing (base 10)
    interSs = logspace(log10(SminmaxnSteps(1)),log10(SminmaxnSteps(2)),NS);
else % use linear spacing
    interSs = linspace(SminmaxnSteps(1),SminmaxnSteps(2),NS);
end

end




