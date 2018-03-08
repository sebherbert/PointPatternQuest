


%% Merge datafiles mother, type 1 and type 2
titidmspotmother.celltypedm = [];

% Add the cell type to the data table
motherType = table(cell(height(titidmspotmother),1),'VariableNames',{'cellType'});
motherType.cellType(:) = {'mother'};
titidmspotmother = horzcat(titidmspotmother,motherType);

type1Type = table(cell(height(titidmspottype1),1),'VariableNames',{'cellType'});
type1Type.cellType(:) = {'type1'};
titidmspottype1 = horzcat(titidmspottype1,type1Type);

type2Type = table(cell(height(titidmspottype2),1),'VariableNames',{'cellType'});
type2Type.cellType(:) = {'type2'};
titidmspottype2 = horzcat(titidmspottype2,type2Type);

% Merge into 1 structure
fullData = vertcat(titidmspotmother,titidmspottype1,titidmspottype2);
% fullData.cellType = char(fullData.cellType);

% Delete empty rows
fullData(isnan(fullData.PositionX(:)),:) = [];

% Update ID for individual ones
fullData.oldID = fullData.ID;
fullData.ID = [1:height(fullData)]';


%% Extract meaningful variables
PARAMS.nbreTp = max(fullData.Time)-(min(fullData.Time)-1);

%% Find duplicates
colors = lines(8);
figure
hold on
dupliMat = [];
for tp = 1:PARAMS.nbreTp
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
    allDistances(allDistances==0) = nan; % duplicates ditances are switched to nan

    % find all the nans in the distance upper matrix
    [fooDupli(:,1), fooDupli(:,2)] = find(isnan(triu(allDistances,1))==1);
    
    % Match cell ID to cell indices
    indices = tpData.ID(fooDupli);
    
    % concatenate with previous timepoints
    dupliMat = vertcat(dupliMat, [fooDupli ones(size(fooDupli,1),1)*tp indices]);    
    

end


%% Display data
colors = lines(8);
figure
hold on
for tp = 1:PARAMS.nbreTp
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


%% Find cell type 
dupliMat = array2table(dupliMat);
dupliMat.Properties.VariableNames = {'cell1','cell2','tp','IDcell1','IDcell2'};

% Yes cell types could be found at previous step but it is less troublesome
% this way...
dupliMat.cellType1 = fullData.cellType(fullData.ID(dupliMat.IDcell1));
dupliMat.cellType2 = fullData.cellType(fullData.ID(dupliMat.IDcell2));

%% Get rid of duplicates
toDeleteList = [];
for bioCell = 1:size(dupliMat,1) % for each duplicate
    if ~(strcmp(dupliMat.cellType1(bioCell),'mother') || strcmp(dupliMat.cellType2(bioCell),'mother'))
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

save('liveDataCurated','fullDataLive')









