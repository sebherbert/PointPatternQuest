

function displayDynNNmap(PARAMS, Rs, Ss, DTfields, statTest, permutPopNames, inNNana)
% all maps of a specific stat test is displayed in a dedicated subplot for
% each deltaT in a 3D surface

figure


numSubP = numSubplots(numel(PARAMS.display.deltaTOIs));

% Display maps for each DeltaTs
for DTfield = PARAMS.display.deltaTOIs(1):PARAMS.display.deltaTOIs(end)

    
    %     displayDynNN(numSubP, DTfields, DTfield, Rs, Ss,...
    %         NNana.perDeltaT.(DTfields{DTfield}).(statTest), statTest);
    
    subplot(numSubP(1),numSubP(2),DTfield)
    hold on
    surf(Rs, Ss, reshape(inNNana.(DTfields{DTfield}).(statTest).map, numel(Ss), numel(Rs)),...
        'EdgeColor', 'None'); %  'FaceColor', 'interp'
    scatter3(inNNana.(DTfields{DTfield}).(statTest).minRval,...
        inNNana.(DTfields{DTfield}).(statTest).minSval,...
        inNNana.(DTfields{DTfield}).(statTest).minVal , 'filled','MarkerFaceColor',[217 83 25]/255);
    title(DTfields{DTfield});
    xlabel('Range'); ylabel('Strength');
    set(gca, 'YScale', 'log');
    legend({sprintf('%s map',statTest) sprintf('Global minimum = %0.3f',...
        inNNana.(DTfields{DTfield}).(statTest).minVal)},...
        'Location', 'northeast');
    colorbar

end

saveas(gcf,sprintf('%smap_permutIn_%s.fig',statTest, permutPopNames));

end







