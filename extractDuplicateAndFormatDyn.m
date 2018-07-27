

function extractDuplicateAndFormatDyn()
% This function merges all type1, type 2 and mother type cells
% it also makes sure that there are no duplicates of each cell in a
% different population. (for example a cell detected as a mother and a type
% 1 at the same time or a cell detected twice)
% It will also add important additional metrics such as a new unique ID and
% reformat the table.


%% Merge datafiles mother, type 1 and type 2

% % load spot files (for example titi_dm_spot_mother.csv)
% % semi-automatic => Changes the type of some columns...
motherFile = uipickfiles('Prompt',... % Mother file
    'Select the data file containing the mother spots (ex: titi_dm_spot_mother.csv) to analyze', 'NumFile', 1);
% Delete empty columns
opts = detectImportOptions(motherFile{1});
spotmother = readtable(motherFile{1},opts);
vars = opts.SelectedVariableNames;
spotmother = spotmother(:,vars);


type1File = uipickfiles('Prompt',... % type1 file
    'Select the data file containing the type1 spots (ex: titi_dm_spot_type1.csv) to analyze', 'NumFile', 1);
% Delete empty columns
opts = detectImportOptions(type1File{1});
spottype1 = readtable(type1File{1},opts);
vars = opts.SelectedVariableNames;
spottype1 = spottype1(:,vars);

type2File = uipickfiles('Prompt',... % type2 file
    'Select the data file containing the type2 spots (ex: titi_dm_spot_type2.csv) to analyze', 'NumFile', 1);
% Delete empty columns
opts = detectImportOptions(type2File{1});
spottype2 = readtable(type2File{1},opts);
vars = opts.SelectedVariableNames;
spottype2 = spottype2(:,vars);

% % load spot files (for example titi_dm_spot_mother.csv)
% % manual => For legacy maintenance => If issue with the format of the
% columns of Unit


spotmother.cellTypeDm = [];

% Add the cell type to the data table
motherType = table(cell(height(spotmother),1),'VariableNames',{'cellType'});
motherType.cellType(:) = {'mother'};
spotmother = horzcat(spotmother,motherType);

type1Type = table(cell(height(spottype1),1),'VariableNames',{'cellType'});
type1Type.cellType(:) = {'type1'};
spottype1 = horzcat(spottype1,type1Type);

type2Type = table(cell(height(spottype2),1),'VariableNames',{'cellType'});
type2Type.cellType(:) = {'type2'};
spottype2 = horzcat(spottype2,type2Type);

% Merge into 1 structure
fullData = vertcat(spotmother,spottype1,spottype2);
% fullData.cellType = char(fullData.cellType);

% Delete empty rows
fullData(isnan(fullData.PositionX(:)),:) = [];

% Update ID for individual ones
fullData.oldID = fullData.ID;
fullData.ID = (1:height(fullData))';


%% Extract meaningful variables
PARAMS.nbreTp = max(fullData.Time)-(min(fullData.Time)-1);

%% Find duplicates
colors = lines(PARAMS.nbreTp);

dupliMat = [];
for tp = min(fullData.Time):max(fullData.Time)
    clear foo*
    tpData = fullData(fullData.Time==tp,:);
    allDistances = squareform(pdist([tpData.PositionX, tpData.PositionY, tpData.PositionZ]));
    
    % => to find the dupliquates of the mother
    %     allDistances = sort(allDistances);
    %     for i = 1:36
    %         temp = fullData.cellType(find(allDistances(:,i)==0));
    %         duplic.(sprintf('tp%d',tp))(i) = temp(2);
    %     end
    
    % => to find all the duplicates
    allDistances(allDistances==0) = nan; % duplicates distances are switched to nan

    % find all the nans in the distance upper matrix
    [fooDupli(:,1), fooDupli(:,2)] = find(isnan(triu(allDistances,1))==1);
    
    % Match cell ID to cell indices
    indices = tpData.ID(fooDupli);
    
    % concatenate with previous timepoints
    dupliMat = vertcat(dupliMat, [fooDupli ones(size(fooDupli,1),1)*tp indices]);    
end


%% Display data
colors = lines(PARAMS.nbreTp);
figure
hold on
for tp = min(fullData.Time):max(fullData.Time)
    tpData = fullData(fullData.Time==tp,:);

    plot3(tpData.PositionX(strcmp(tpData.cellType,'mother')),...
        tpData.PositionY(strcmp(tpData.cellType,'mother')),...
        tpData.PositionZ(strcmp(tpData.cellType,'mother')),...
        '^','Color',colors(tp,:));
    plot3(tpData.PositionX(strcmp(tpData.cellType,'type1')),...
        tpData.PositionY(strcmp(tpData.cellType,'type1')),...
        tpData.PositionZ(strcmp(tpData.cellType,'type1')),...
        '.','Color',colors(tp,:));
    plot3(tpData.PositionX(strcmp(tpData.cellType,'type2')),...
        tpData.PositionY(strcmp(tpData.cellType,'type2')),...
        tpData.PositionZ(strcmp(tpData.cellType,'type2')),...
        'o','Color',colors(tp,:));
end


%% Find cell type of the duplicates
dupliMat = array2table(dupliMat);
dupliMat.Properties.VariableNames = {'cell1','cell2','tp','IDcell1','IDcell2'};

% Yes cell types could be found at previous step but it is less troublesome
% this way...
dupliMat.cellType1 = fullData.cellType(fullData.ID(dupliMat.IDcell1));
dupliMat.cellType2 = fullData.cellType(fullData.ID(dupliMat.IDcell2));

%% Get rid of duplicates
toDeleteList = [];
for bioCell = 1:size(dupliMat,1) % for each duplicate
    if ~((strcmp(dupliMat.cellType1(bioCell),'mother') && strcmp(dupliMat.cellType2(bioCell),'type2')) || ...
            (strcmp(dupliMat.cellType1(bioCell),'type2') && strcmp(dupliMat.cellType2(bioCell),'mother')))
        if strcmp(dupliMat.cellType1(bioCell),dupliMat.cellType2(bioCell))
            toDeleteList = [toDeleteList dupliMat.IDcell1(bioCell)];
        else
            fprintf('Warning! cells IDs %d and %d are dupliquates of type %s and %s! Please pick one.\n',...
                dupliMat.IDcell1(bioCell), dupliMat.IDcell2(bioCell), ...
                dupliMat.cellType1{bioCell},dupliMat.cellType2{bioCell});
            return
        end
    end
end
fprintf('Deleting cells ID %d\n', toDeleteList);
fullData(fullData.ID(toDeleteList),:) = [];

fullDataLive = fullData;

fullDataLive(:,{'Category','Collection'}) = [];

save('liveDataCurated','fullDataLive');

end







