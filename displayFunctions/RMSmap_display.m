
%% explore RMS map version v2p1

% load ...dataCombined_RMSMap
% CAUTION, will load the whole file into MATLAB, easily eating 60Go of RAM

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
            dataCombined_RMSMap.(fieldList{field}).t3vst2.medRMS));
    end
end
variableNames = {'Range','Strength','finalRMS'};
allModels.Properties.VariableNames = variableNames;

% reshape final RMSE into a single 2D matrix
RMSEmap = reshape(allModels.finalRMS,[numel(unique(allModels.Strength)),numel(unique(allModels.Range))]);
Rangemap = reshape(allModels.Range,[numel(unique(allModels.Strength)),numel(unique(allModels.Range))]);
Strengthmap = reshape(allModels.Strength,[numel(unique(allModels.Strength)),numel(unique(allModels.Range))]);


% surf(Rs, Ss, statMap, 'EdgeColor', 'None', 'FaceColor', 'interp');
surf(Rangemap, Strengthmap, RMSEmap, 'EdgeColor', 'None', 'FaceColor', 'interp');

xlabel('Range (\mum)'); ylabel('Strength');
set(gca, 'YScale', 'log');

colorbar

az = 0;
el = 90;
view(az, el);



