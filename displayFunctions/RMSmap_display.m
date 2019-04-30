

function RMSmap_display(dataCombined_RMSMap, saveFilePath)
%% explore RMS map version v2p1

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
            dataCombined_RMSMap.(fieldList{field}).t3vst3.medRMS));
    end
end
variableNames = {'Range','Strength','finalRMS'};
allModels.Properties.VariableNames = variableNames;

% reshape final RMSE into a single 2D matrix
RMSEmap = reshape(allModels.finalRMS,[numel(unique(allModels.Strength)),numel(unique(allModels.Range))]);
Rangemap = reshape(allModels.Range,[numel(unique(allModels.Strength)),numel(unique(allModels.Range))]);
Strengthmap = reshape(allModels.Strength,[numel(unique(allModels.Strength)),numel(unique(allModels.Range))]);


% surf(Rs, Ss, statMap, 'EdgeColor', 'None', 'FaceColor', 'interp');
figure
surf(Rangemap, Strengthmap, RMSEmap, 'EdgeColor', 'None', 'FaceColor', 'interp');

xlabel('Range (\mum)'); ylabel('Strength');
set(gca, 'YScale', 'log');

colorbar

az = 0;
el = 90;
view(az, el);

saveas(gcf,sprintf('%s_t3vst3_RMSEmap', saveFilePath));
saveas(gcf,sprintf('%s_t3vst3_RMSEmap.png', saveFilePath));

end



