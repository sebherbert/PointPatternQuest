

function extract_data_bioCompare()
% This function will fetch the meaning values in all of the possible
% samples and merge them together for simpler metric visualisation



fetchError = 1; % Use if error estimation has already been done


%% load all data
% Data must follow the .mat file example of
% 170126_2month_20x_subdiv_L_2nodb_noDl_xyzCASE1_allSample_fittedModel.mat
fileToOpen = uipickfiles('Prompt','Please, select the correct file to analyse','FilterSpec','*.mat','REFilter','fit','REDirs',0);
% temp = load('loader.mat');
% fileToOpen = temp.fileToOpen;
fullTable = table;
if fetchError
    errorTable = table;
end

rootFolder = [];

% statTests = {'t2vst2','t3vst3','t3vst2'};
statTests = {'t3vst2'};

for file = 1:length(fileToOpen)
    load(fileToOpen{file})
    [path, fileName, ~] = fileparts(fileToOpen{file});
    
    if ~strcmp(rootFolder,path)
        folderIndexes = strfind(path,filesep) ;
        lastFolderName = path(folderIndexes(end)+1:end);
        fprintf('\nNew folder: %s\n',lastFolderName);
        rootFolder = path;
    end
    
    if contains(path,'2to3month')
        bioCondition = 'Young';
    elseif contains(path,'4to3month')
        bioCondition = 'Old';
    elseif contains(path,'LY24h')
        bioCondition = 'Drug';
    else
        fprintf('ERROR in biological condition parsing\n');
        return
    end
        
    if contains(fileName,'_allSample_')
        brainPart = 'allSample';
    elseif contains(fileName,'_dm_')
        brainPart = 'dm';
    elseif contains(fileName,'_da_')
        brainPart = 'da';
    elseif contains(fileName,'_dl_')
        brainPart = 'dl';
    else
        fprintf('ERROR in brain part parsing\n');
        return
    end
        
    fprintf('file %s\toptiR0 = %d\toptiS0 = %d\n',lastFolderName,...
        dataCombined.(statTests{1}).PARAMS.optiR0,...
        dataCombined.(statTests{1}).PARAMS.optiS0);
    
    for statTest = 1 : numel(statTests)
        fullTable = vertcat(fullTable,cell2table({lastFolderName, ...
            bioCondition, brainPart, statTests{statTest}...
            dataCombined.(statTests{statTest}).fit.Range, ...
            dataCombined.(statTests{statTest}).fit.Strength, ...
            dataCombined.(statTests{statTest}).fit.medRMS, ...
            dataCombined.(statTests{statTest}).PARAMS.numPermut, ...
            dataCombined.(statTests{statTest}).PARAMS.fitMaxIter, ...
            dataCombined.(statTests{statTest}).PARAMS.optiR0, ...
            dataCombined.(statTests{statTest}).PARAMS.optiS0}));
    end
    
    if fetchError
        for statTest = 1 : numel(statTests)
            errorTable = vertcat(errorTable,cell2table({real(dataCombined.(statTests{statTest}).assoDer.fast.RMS_Rferror), ...
                dataCombined.(statTests{statTest}).assoDer.fast.RMS_Sferror, ...
                dataCombined.(statTests{statTest}).assoDer.fast.RMS2nd_Rf, ...
                dataCombined.(statTests{statTest}).assoDer.fast.RMS2nd_Sf}));
            if errorTable.Var3(end) < 0
                errorTable.Var1(end) = nan;
            end
        end
    end
end

variableNames = {'Sample','bioCondition','brainPart','distTest','Range','Strength','finalRMS','numPermut','fitMaxIter','R0','S0'};
fullTable.Properties.VariableNames = variableNames;

if fetchError
    errorTable.Properties.VariableNames = {'RMS_Rferror', 'RMS_Sferror', 'RMSf2nd_RR', 'RMSf2nd_SS'};
    fullTable = horzcat(fullTable,errorTable);
end


%% Parse and plot data in a subplot

fieldsToPlot = {'Range','Strength','finalRMS'};
lineAxis = 'distTest';
% allBrainParts = {'allSample','da','dm'};
allBrainParts = {'dm'};
allBioConditions = {'Young','Old','Drug'};
scatterDotStyle = {'v','o','square'};

figure
subPlotNbr = 0;
for lineValNth = 1:numel(statTests)
    % find the test of interest
    lineVal = statTests{lineValNth};
    
    % extract the test of interest
    dataForPlotInLine = fullTable(strcmp(fullTable.(lineAxis),lineVal),:);

    % for each field
    for field = 1:numel(fieldsToPlot)
        fdName = fieldsToPlot{field};
        subPlotNbr = subPlotNbr + 1;
        subplot(numel(statTests),numel(fieldsToPlot),subPlotNbr)
        hold on
        xPos = 0;
        
        % for each brainPart
        for bp = 1:numel(allBrainParts)
            bpName = allBrainParts{bp};
            
            % for each bioCondition
            for bc = 1:numel(allBioConditions)
                bcName = allBioConditions{bc};
                
                xPos = xPos + 1;
                scatterDot = scatterDotStyle{bc};
                
                % parse full data for data of interest
                dataForSubPlot = dataForPlotInLine...
                    (strcmp(dataForPlotInLine.brainPart,bpName),:);
                dataForSubPlot = dataForSubPlot...
                    (strcmp(dataForSubPlot.bioCondition,bcName),:);
                
                % Check if data exist for this specific set of parameters
                if numel(dataForSubPlot)==0
                    disp('no data');
                    continue
                end
                
                % display
                h(bc) = scatter(xPos.*ones(numel(dataForSubPlot.Range),1),...
                    dataForSubPlot.(fdName),50,...
                    dataForSubPlot.finalRMS,scatterDot,'LineWidth', 2);
            end
            xPos = xPos+1;
        end
        
        % Improve subplot display
        title(sprintf('%s (%s)',fdName,lineVal));

        xticks(2:4:12);
        xticklabels(allBrainParts)

        current = gca;
        
        if strcmp(fdName,'Strength')
            set(current, 'YScale', 'log')
            current.YAxis.Limits = [0.1 50];

        elseif strcmp(fdName,'finalRMS')
            legend(h,allBioConditions,'Location','eastoutside');
            legend boxoff;
        
            c = colorbar;
            c.Label.String = 'RMS';
            
            current.XAxis.Limits = [0 12];
        end
        grid on
    end
end

end




