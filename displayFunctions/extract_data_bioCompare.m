

function extract_data_bioCompare

fetchError = 0; % Use if error estimation has already been done


%% load all data
% Data must follow the .mat file example of
% 170126_2month_20x_subdiv_L_2nodb_noDl_xyzCASE1_allSample_fittedModel.mat
fileToOpen = uipickfiles('Prompt',['Please, select all the correct files to analyse (folder name must'...
    ' be the condition name)'],'FilterSpec','*.mat','REFilter','fit','REDirs',0);
% temp = load('loader.mat');
% fileToOpen = temp.fileToOpen;
fullTable = table;
if fetchError
    errorTable = table;
end

rootFolder = [];

statTests = {'t2vst2','t3vst3','t3vst2'};
% statTests = {'t3vst2'};

for file = 1:length(fileToOpen)
    load(fileToOpen{file})
    [path, fileName, ~] = fileparts(fileToOpen{file});
    data1File = dataCombined_fitModel.data;
    
    fprintf('\nLoading file %s\n', fileName);
    sampleName = regexprep(fileName,'_xyzCASE.*','');
    
    %     if ~strcmp(rootFolder,path)
    %         folderIndexes = strfind(path,filesep) ;
    %         sampleName = path(folderIndexes(end)+1:end);
    %         fprintf('\nNew folder: %s\n',sampleName);
    %         rootFolder = path;
    %     end
    
    % find the experimental condition
    bioCondition = regexprep(path,'.*/','');
    fprintf('bioCondition = %s\n', bioCondition);
    
    
    %  => old version with fixed names to find in file name
    %     if contains(path,'2to3month')
    %         bioCondition = 'Young';
    %     elseif contains(path,'4to3month')
    %         bioCondition = 'Old';
    %     elseif contains(path,'LY24h')
    %         bioCondition = 'Drug';
    %     else
    %         fprintf('ERROR in biological condition parsing\n');
    %         return
    %     end
        
    % find the brainpart
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
    fprintf('brainPart = %s\n', brainPart);
        
    fprintf('optiR0 = %d\toptiS0 = %d\n',...
        data1File.(statTests{1}).PARAMS.optiR0,...
        data1File.(statTests{1}).PARAMS.optiS0);
    
    for statTest = 1 : numel(statTests)
        fullTable = vertcat(fullTable,cell2table({sampleName, ...
            bioCondition, brainPart, statTests{statTest}...
            data1File.(statTests{statTest}).fit.Range, ...
            data1File.(statTests{statTest}).fit.Strength, ...
            data1File.(statTests{statTest}).fit.medRMS, ...
            data1File.(statTests{statTest}).PARAMS.numPermut, ...
            data1File.(statTests{statTest}).PARAMS.fitMaxIter, ...
            data1File.(statTests{statTest}).PARAMS.optiR0, ...
            data1File.(statTests{statTest}).PARAMS.optiS0}));
    end
    
    if fetchError
        for statTest = 1 : numel(statTests)
            errorTable = vertcat(errorTable,cell2table({real(data1File.(statTests{statTest}).assoDer.fast.RMS_Rferror), ...
                data1File.(statTests{statTest}).assoDer.fast.RMS_Sferror, ...
                data1File.(statTests{statTest}).assoDer.fast.RMS2nd_Rf, ...
                data1File.(statTests{statTest}).assoDer.fast.RMS2nd_Sf}));
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
allBrainParts = unique(fullTable.brainPart);
allBioConditions = unique(fullTable.bioCondition, 'stable');
scatterDotStyle = {'v','o','s','^','d','+'};

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

        current = gca;

        
        if numel(allBrainParts)==1
            xticks(1:numel(allBioConditions))
            xticklabels(allBioConditions);
            xtickangle(45);
            current.XAxis.Limits = [0 numel(allBioConditions)+1];
        else
            xticks(2:4:numel(allBrainParts)*4); % 4 could be swaped with number of biocondition+1 instead?
            xticklabels(allBrainParts);
        end
        
        if strcmp(fdName,'Strength')
            set(current, 'YScale', 'log')
            current.YAxis.Limits = [0.1 50];

        elseif strcmp(fdName,'finalRMS')
            legend(h,allBioConditions,'Location','eastoutside');
            legend boxoff;
        
            c = colorbar;
            c.Label.String = 'RMS';
            if numel(allBrainParts)==1
                current.XAxis.Limits = [0 numel(allBioConditions)+1];
            else
                current.XAxis.Limits = [0 12]; % idem higher
            end
        end
        grid on
    end
end

end



