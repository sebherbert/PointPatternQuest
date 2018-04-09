

function RMSmap_display()
% DEPRECATED: This function displays the RMS map of any matrix of R and S values as a
% 2D plot.
% version v1p1 => Currently not working with variable input size... 

% Select a file containing th RMSE map (named dataCombinedModels)
fileToOpen = uipickfiles('Prompt','Please, select the correct file to analyse','FilterSpec','*.mat','REDirs',0,'numfiles',1);

% Load the RMSE map file
load(fileToOpen{1});

% get list of models
fieldList = fieldnames(dataCombinedModels);

% get table of results
allModels = table;
for field = 1:numel(fieldList) % for each field
    % Check if field is a model
    if contains(fieldList{field},'Model')
        allModels = vertcat(allModels, ...
            table(dataCombinedModels.(fieldList{field}).effectRange, ...
            dataCombinedModels.(fieldList{field}).effectStrength, ...
            dataCombinedModels.(fieldList{field}).t3vst2.medRMS));
    end
end
variableNames = {'Range','Strength','finalRMS'};
allModels.Properties.VariableNames = variableNames; 


%% Display RMS map
figure
X = [1,1];
Y = [numel(unique(allModels.Strength)), numel(unique(allModels.Range))];

imagesc(X, Y, allModels.finalRMS)

map = colormap(parula(1024));
minmaxRMS = (allModels.finalRMS-min(allModels.finalRMS))/(max(allModels.finalRMS)-min(allModels.finalRMS));
imshow(reshape(minmaxRMS,[Y(1),Y(2)])*length(map),map)

scatter(allModels.Range,allModels.Strength,90,allModels.finalRMS,'filled');

c = colorbar;
c.Label.String = 'RMS';
caxis([0.02 0.5])

hi = imagesc( c )
xt = get(gca, 'XTick' )
xs = num2str(x(xt))
set(gca, 'XTickLabel', xs)

yt = get(gca, 'YTick')
ys = num2str(y(yt))
set(gca, 'YTickLabel', ys)

set( gca, 'XTick', get(gca, 'XTick') + 1)
set( gca, 'XTickLabel', num2str( x( get(gca, 'XTick') ) ) )


end

