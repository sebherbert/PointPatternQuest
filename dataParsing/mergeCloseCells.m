%%
% 
%  This function allow the user to merge pairs of cells are too close from
%  each others. It uses the NNExp structure and the parameters
%  - PARAMS.cellMerge.cellType for the cell types to be merged and
%  - PARAMS.cellMerge.distThresh for the maximum distance below which objects
%  should be merged.
% 


function NNExp = mergeCloseCells(PARAMS, NNExp)

for popOIs = 1:length(PARAMS.cellMerge.cellType)
    cellType2merge = PARAMS.cellMerge.cellType{popOIs}; % specify cell type
    % Simplify whole matrix to only keep the cell type of interest
    cells2check = NNExp(NNExp.cellType == cellType2merge,:);
    % Keep in check the max ID for additional merged cells
    maxID = max(NNExp.cellID);
    
    % Measure all the distances between the cells
    dists = pdist(cells2check.pos3D, 'euclidean');
    
    % Recreate the square form of the linear indexes and keep only a single occurence of each minimum
    distSqBL = tril(squareform(dists)); 

    % Only keep indexes for which dist<threshold
    distSqBL(distSqBL>PARAMS.cellMerge.distThresh) = 0;
    [row, col, values] = find(distSqBL);
    
    % Sort positions by smallest distances
    [sortedValues, idxValues] = sort(values);
    row = row(idxValues);
    col = col(idxValues);
    values = values(idxValues);
    
    % Find cell's ID associated to each pair
    cellPairsID = [cells2check.cellID(row), cells2check.cellID(col)];
    
    % for each pair
    mergedCells = cells2check(1,:);
    mergedCells(1,:) = [];
    
    % Create new table to dump the new merged cells
    mergedCells = cell2table(cell(0,width(NNExp)));
    mergedCells.Properties.VariableNames = NNExp.Properties.VariableNames;

    % Create new table to dump the deleted cells
    deletedCells = cell2table(cell(0,width(NNExp)));
    deletedCells.Properties.VariableNames = NNExp.Properties.VariableNames;
    
    for pair = 1:size(cellPairsID,1)
        presentPair = cellPairsID(pair, :);
        % check if both ID are still present
        if (isempty(cells2check(cells2check.cellID == presentPair(1),:)) || ...
                isempty(cells2check(cells2check.cellID == presentPair(2),:)))
            continue % at least one has already been merged => Continue
        else
            % Add new merged cell
            pos3D_1 = cells2check.pos3D(cells2check.cellID == presentPair(1), :);
            pos3D_2 = cells2check.pos3D(cells2check.cellID == presentPair(2), :);
            %             pdist2(pos3D_1, pos3D_2)

            % Prepare new cell table holder
            newPos3D = cell2table(cell(1,width(mergedCells)));
            newPos3D.Properties.VariableNames = mergedCells.Properties.VariableNames;
            
            % Fill up new cell table holder
            newPos3D.cellID = maxID + height(mergedCells) + 1;
            newPos3D.cellType = cellType2merge;
            newPos3D.nearestNeighbour = nan;
            newPos3D.pos3D = mean([pos3D_1 ; pos3D_2]);
            
            % Merge it with rest of new cells
            mergedCells = [mergedCells ; newPos3D];
            
            % Copy lines  from original table into a temporary file
            deletedCells = [deletedCells ; cells2check(cells2check.cellID == presentPair(1),:)];
            deletedCells = [deletedCells ; cells2check(cells2check.cellID == presentPair(2),:)];
            
            % Copy lines into a temporary file and delete them from original table
            cells2check(cells2check.cellID == presentPair(1),:) = [];
            cells2check(cells2check.cellID == presentPair(2),:) = [];
        end
    end
    
    % Merge back into population of interest
    cells2check = [cells2check ; mergedCells];
    
    % Merge back into overall population
    NNExp(NNExp.cellType == cellType2merge,:) = [];
    NNExp = [NNExp; cells2check];
    NNExp.cellID = [1:height(NNExp)]';
    
end

% Update nearest neighbour field of the table
newNNmat = squareform(pdist(NNExp.pos3D));
newNNmat(newNNmat==0) = nan;
NNExp.nearestNeighbour = min(newNNmat)';

end
