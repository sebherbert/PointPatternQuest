
function displayModelSerie(dataCombinedModels, PARAMS)
% Display all overlayed models with experimental data as cdfs

% List tested models
fooFields = fieldnames(dataCombinedModels);
models = {};
for field = 1:numel(fooFields)
    if contains(fooFields{field},'Model_')
        models{numel(models)+1} = fooFields{field};
    end
end

% For each statistical test => Assumes all models are tested for the same
% statistical tests in terms crossed population
fooDispTests = fieldnames(dataCombinedModels.(models{1}));
dispTests = {};
for field = 1:numel(fooDispTests)
    if regexp(fooDispTests{field},'t[0-9]vst[0-9]')
        dispTests{numel(dispTests)+1} = fooDispTests{field};
    end
end

clear foo*

% For each dispersion test found (usually 3)
for dispTest = 1:numel(dispTests)
    clear h legs
    h = [];
    legs = {};
    
    figure
    hold on

    % Prepare figure colors for each model
    %     colors = [jet(numel(models)+1) ones(numel(models)+1,1)*0.25];
    
    % use first model to extract the experimental data
    currentData = dataCombinedModels.(models{1}).(dispTests{dispTest});
    colors = [lines(2) ones(2,1)*0.25]; % => Keep the 2nd default color for
    % the experimental dataset   
    h(length(h)+1) = displayIndivCDF(currentData.expCDFs.x,currentData.expCDFs.f,...
        colors(2,:),currentData.expCDFs.f5,currentData.expCDFs.f95,'.-');
    legs{1} = 'Experimental data';
    % Prepare figure colors for each model
    %     colors(2,:) = [];
    colors = [winter(numel(models)) ones(numel(models),1)*0.25];
    % Loop the display of the modeled and experimental CDFs
    for model = 1:numel(models)
        currentData = dataCombinedModels.(models{model}).(dispTests{dispTest});
        
        h(length(h)+1) = displayIndivCDF(currentData.simuCDFs.x,currentData.simuCDFs.f50pc,...
            colors(model,:),currentData.simuCDFs.f5pc,currentData.simuCDFs.f95pc,'-');
        legs{numel(legs)+1} = regexprep(models{model},'_',' ');        
    end
    
    % Massage display
    ylabel('Cumulative cell frequency');
    xlabel('Distance to nearest neighbor (µm)');
    title(regexprep(sprintf('%s - %s - %s', PARAMS.name, PARAMS.brainPart, dispTests{dispTest}),'_',' '));
    legend(h,legs,'Location','SouthEast');
    PARAMS.axis = [0 100 0 1];

    saveas(gcf,regexprep(sprintf('%s - %s - %s - All models', PARAMS.name, PARAMS.brainPart, dispTests{dispTest}),'_',' '));
    
end


end

function h = displayIndivCDF(x,f,color,f5,f95,lineType)

h = plot(x,f,lineType,'linewidth',1.5,'Color',color(1,1:3));
plot(x,f5,'--','linewidth',1,'Color',color);
plot(x,f95,'--','linewidth',1,'Color',color);

end