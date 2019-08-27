

function RMSmap_display(dataCombined_RMSMap, saveFilePath, popTestsFields)
%% explore RMS map version v2p1


for popTestsNum = 1:length(popTestsFields)
    
    popTests = popTestsFields{popTestsNum};
    
    % get list of models
    fieldList = fieldnames(dataCombined_RMSMap);
    
    % get table of results
    allModels = table;
    for field = 1:numel(fieldList) % for each field
        % Check if field is a model
        if contains(fieldList{field},'Model')
            allModels = vertcat(allModels, ...
                table(dataCombined_RMSMap.(fieldList{field}).effectRange, ...
                dataCombined_RMSMap.(fieldList{field}).effectStrength, ...
                dataCombined_RMSMap.(fieldList{field}).(popTests).medRMS));
        end
    end
    variableNames = {'Range','Strength','finalRMS'};
    allModels.Properties.VariableNames = variableNames;
    
    % reshape final RMSE into a single 2D matrix
    RMSEmap = reshape(allModels.finalRMS,[numel(unique(allModels.Strength)),numel(unique(allModels.Range))]);
    Rangemap = reshape(allModels.Range,[numel(unique(allModels.Strength)),numel(unique(allModels.Range))]);
    Strengthmap = reshape(allModels.Strength,[numel(unique(allModels.Strength)),numel(unique(allModels.Range))]);
    
    saveMapImage(Rangemap, Strengthmap, RMSEmap, saveFilePath, popTests)
    
    % saveMapTable(Rangemap, Strengthmap, RMSEmap, saveFilePath, popTests)
end

end


function saveMapImage(Rangemap, Strengthmap, RMSEmap, saveFilePath, popTests)
% Display and save the map image
% surf(Rs, Ss, statMap, 'EdgeColor', 'None', 'FaceColor', 'interp');
figure
surf(Rangemap, Strengthmap, RMSEmap, 'EdgeColor', 'None', 'FaceColor', 'interp');

xlabel('Range (\mum)'); ylabel('Strength');
set(gca, 'YScale', 'log');

colorbar

az = 0;
el = 90;
view(az, el);

saveas(gcf,sprintf('%s_%s_RMSEmap', saveFilePath, popTests));
saveas(gcf,sprintf('%s_%s_RMSEmap.png', saveFilePath, popTests));

end


function saveMapTable(Rangemap, Strengthmap, RMSEmap, saveFilePath, popTests)
% save the map table

% reformat the RMSEmap as a clean table
mapTable = table();



end
