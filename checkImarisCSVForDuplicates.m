

function checkImarisCSVForDuplicates()
% Check that the CSVs produced by Imaris are not displaying any duplicate


% select the different types for the different parts
fileToOpen = uipickfiles('Prompt','Select the csv files of interest (types1,2,3) (da and dm)');

filePaths = {};
for files = 1:numel(fileToOpen)
    [filepath,name,ext] = fileparts(fileToOpen{files});
    filePaths{files}.filepath = filepath;
    filePaths{files}.name = name;
    filePaths{files}.ext = ext;
end
clear filepath name ext

% Structure the data
parts = {'da' 'dm'};
cellTypes = {'type1' 'type2' 'type3' 'type3_notSOX2'};
data = {};
tableLegs = {'xPos' 'yPos' 'zPos' 'time'};
for files = 1:numel(fileToOpen)
    tempTable = {};
    tempTable = readtable(fileToOpen{files});
    if isempty(regexpi(filePaths{files}.name,'not'))
        fprintf('Parsing: %s\n',filePaths{files}.name);
        cellType = [];
        % Checking each populations combined
        for part = 1:numel(parts)
            if regexp(filePaths{files}.name,parts{part})
                break
            end
        end
        for cellType = 1:numel(cellTypes)
            if regexp(filePaths{files}.name,cellTypes{cellType})
                break
            end
        end
        fprintf('is file %s , %s\n',cellTypes{cellType},parts{part})
        data.(parts{part}).(cellTypes{cellType}) = table(tempTable.PositionX, tempTable.PositionY, tempTable.PositionZ,...
            tempTable.Time, 'VariableNames',tableLegs);
    else
        %         fprintf('NOT parsing: %s\n',filePaths{files}.name);
        if regexp(filePaths{files}.name,'type3')
            fprintf('Parsing: %s\n',filePaths{files}.name);
            cellType = 'type3_notSOX2';
            if regexp(filePaths{files}.name,'_da')
                data.da.(cellType) = table(tempTable.PositionX, tempTable.PositionY, tempTable.PositionZ, tempTable.Time, ...
                    'VariableNames',tableLegs);
                fprintf('is file %s , %s\n',cellType,'da');
            elseif regexp(filePaths{files}.name,'_dm')
                data.dm.(cellType) = table(tempTable.PositionX, tempTable.PositionY, tempTable.PositionZ, tempTable.Time, ...
                    'VariableNames',tableLegs);
                fprintf('is file %s , %s\n',cellType,'dm');
            end
        else
            fprintf('WHAAAAAT???... Call Sebastien, the REGEXP are inadapted\n')
        end
    end
end



% Check for duplicates
parts = {'da' 'dm'};
cellTypes = {'type1' 'type2' 'type3' 'type3_notSOX2'};
% Checking each populations combined
for part = 1:numel(parts)
    for cellType1 = 1:numel(cellTypes)
        for cellType2 = cellType1+1:numel(cellTypes)
            full = [data.(parts{part}).(cellTypes{cellType1}); data.(parts{part}).(cellTypes{cellType2})];
            noDuplicate = unique(full);
            if numel(full) == numel(noDuplicate)
                fprintf('There are no duplicates between %s and %s in %s\n', cellTypes{cellType1}, ...
                    cellTypes{cellType2}, parts{part});
            else
                fprintf('There are %d duplicates between %s and %s in %s\n', ...
                    size(full,1)-size(noDuplicate,1), cellTypes{cellType1}, cellTypes{cellType2}, parts{part});
            end
        end
    end
end

display3Dpositions(data, parts, cellTypes)
    
end

function display3Dpositions(data, parts, cellTypes)
% this function is coded for up to 4 types in n parts


colTable = [0 113 255 ; 191 0 191 ; 236 176 31 ; 0 127 0]/255;
for part = 1:numel(parts)
    figure
    hold on
    title(sprintf('part %s of the brain',parts{part}));
    for cellType = 1:numel(cellTypes)
        dataToPlot = data.(parts{part}).(cellTypes{cellType});
        plot3(dataToPlot.xPos , dataToPlot.yPos , dataToPlot.zPos,'.','Color',colTable(cellType,:));
    end
    legend(cellTypes)
    axis equal
end


end








