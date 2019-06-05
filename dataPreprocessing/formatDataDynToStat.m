%{
This function takes a fullDataLive input and outputs a new set of data
files that are fit for a static analysis

one new folder for each sample 
each folder contains a separated .mat file for each timepoint

required variables in the new .mat file
- title must contain '_xyzCASE1' : to specify the brain part
- S: Cell type number (1xN double)
- d123_1: 3D euclidian dist to closest neighbour (1xN double) in the whole population
- x, y, z : x ,y ,z positions of the nuclei (in pix) (Nx1 single)

%}




function formatDataDynToStat()

fileToOpen = uipickfiles('Prompt',...
    'Select the correct file(s) to process (example: liveDataCurated.mat)');

% new folder name => dump folder
newFold = 'static';

for fileOI=1:length(fileToOpen)
    format1file(fileToOpen{fileOI}, newFold);    
end

end


function format1file(fullFilePath, newFold)

[filePath, fileNameOri, ~]  = fileparts(fullFilePath);

liveShape = load(fullFilePath);

filePath = strcat(filePath, filesep, newFold);

mkdir(filePath);



% for each time point
for tp = 1:max(liveShape.fullDataLive.Time)-min(liveShape.fullDataLive.Time)+1
    temptp = liveShape.fullDataLive(liveShape.fullDataLive.Time == tp, :);
    
    [S, d123_1, x, y, z] = format1tp(temptp);
       
    % save tp and clear it
    brainPart = '_xyzCASE1';
    fileNameNew = strcat(filePath, filesep, fileNameOri, brainPart, '_tp', num2str(tp));
    save(fileNameNew , 'S', 'd123_1', 'x', 'y', 'z');
    clear formattedTp
end

end



function [S, d123_1, x, y, z] = format1tp(temptp)

% % Initialize output struct
% formattedTp = {};

% All mother cells are duplicated in type2 => get rid of them
temptp(strcmp(temptp.cellType,'mother'),:) = [];

% Initialize output struct
S = 1:height(temptp);
d123_1 = 1:height(temptp);

for bioCell = 1:height(temptp)
    % Cell type
    if strcmp(temptp.cellType(bioCell),'mother') 
        S(bioCell) = 0; % => should be none since all deleted
    elseif strcmp(temptp.cellType(bioCell),'type1') 
        S(bioCell) = 1;
    elseif strcmp(temptp.cellType(bioCell),'type2')
        S(bioCell) = 2;
    else
        error('Unrecognized cellType of cell ID %i\nExiting\n', temptp.ID(bioCell));
    end
end

% x,y,z position
x = single(temptp.PositionX);
y = single(temptp.PositionY);
z = single(temptp.PositionZ);

% 3D dist
d123_1 = pdist2( [x,y,z], [x,y,z], 'euclidean', 'Smallest', 2);
d123_1 = d123_1(2,:);

% % format output
% formattedTp.S = S;
% formattedTp.d123_1 = d123_1;
% formattedTp.x = x;
% formattedTp.y = y;
% formattedTp.z = z;

end





% total = 0;
% for tp = 1:8
%     NbrI = 0;
%     for i = 1:height(liveShape.fullDataLive)
%         if liveShape.fullDataLive.Time(i)==tp
%             NbrI = NbrI+1;
%         end
%     end
%     fprintf(' at tp=%i : NbrI=%i\n',tp,NbrI);
%     total = NbrI+total;
% end
    