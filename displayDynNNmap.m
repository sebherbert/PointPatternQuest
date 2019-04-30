

function displayDynNNmap(PARAMS, Rs, Ss, DTfields, statTest, permutPopNames, inNNana, doIso)
% all maps of a specific stat test is displayed in a dedicated subplot for
% each deltaT in a 3D surface or an isosurface

figure

numSubP = numSubplots(numel(PARAMS.display.deltaTOIs));

if doIso
    isoStatus = 'iso';
else
    isoStatus = '';
end

% Display maps for each DeltaTs
for DTfield = PARAMS.display.deltaTOIs(1):PARAMS.display.deltaTOIs(end)
    
    statMap = reshape(inNNana.(DTfields{DTfield}).(statTest).map, numel(Ss), numel(Rs));
    
    subplot(numSubP(1),numSubP(2),DTfield)
    hold on
    
    scatter3(inNNana.(DTfields{DTfield}).(statTest).minRval,...
        inNNana.(DTfields{DTfield}).(statTest).minSval,...
        inNNana.(DTfields{DTfield}).(statTest).minVal , 'filled','MarkerFaceColor',[217 83 25]/255);
    
    if doIso
        contour(Rs, Ss, statMap);
    else
        %         surf(Rs, Ss, statMap, 'EdgeColor', 'None'); %  'FaceColor', 'interp'
        surf(Rs, Ss, statMap, 'EdgeColor', 'None', 'FaceColor', 'interp'); 
    end
    
    legend({sprintf('GlobMin= %0.3f', inNNana.(DTfields{DTfield}).(statTest).minVal)...
        sprintf('%s %smap',statTest, isoStatus)}, 'Location', 'southoutside');

    title(sprintf('%s, N=%d',DTfields{DTfield},inNNana.(DTfields{DTfield}).exp.Ncells));
    xlabel('Range (\mum)'); ylabel('Strength');
    set(gca, 'YScale', 'log');

    colorbar

end

if ~PARAMS.dummy
    saveas(gcf,sprintf('%s%smap_permutIn_%s.fig',statTest, isoStatus, permutPopNames));
end

end







