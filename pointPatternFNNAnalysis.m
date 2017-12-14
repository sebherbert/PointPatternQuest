


function fullResults = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS)
% Measure and compare the nearest neighbour distance in experimental data
% and simulations
%{
Inputs:
  fullPath = absolute path plus new file name root
  NNExp = table of : cellType, NN distance, 3D positions of every cell
  cCDists = all cell to cell distances (square array of NxN cells)
  pops.popSource = celltype of the pop source of a possible effect
  pops.popTarget = celltype of the pop affected by the pop source
  pops.popPermut = celltypes among which the random permutations are done
  PARAMS = set of parameters (such as binSize or number of simulations)
Output:
  fullResults = structure containing the NN distribution for the simulations
  as well as the experimental data

%}

% % For simpler reading
% cellType = NNExp.cellType;
% r = PARAMS.binSize;

% Find the average cell size by averaging the minimum cell to cell distance
CellDiameter = mean(NNExp.nearestNeighbour);

%% Extract the 3D positions of the 2 cell populations of interest
popSource3Dpos = table2array(NNExp(NNExp.cellType == pops.popSource, {'pos3D'}));

popTarget3Dpos = table2array(NNExp(NNExp.cellType == pops.popTarget, {'pos3D'}));
nTarget = length(popTarget3Dpos);

rowPermut = logical(sum(NNExp.cellType == pops.popPermut,2));
popPermut3Dpos = table2array(NNExp(rowPermut, {'pos3D'}));

if pops.popTarget == pops.popSource
    samePop = true;
else 
    samePop = false;
end

% Find nearest neighbour
dnExp = findNN(popSource3Dpos, popTarget3Dpos, samePop, pops)';

% Display experimental populations cells
% figure
% hold on
% plot3(popPermut3Dpos(:,1),popPermut3Dpos(:,2),popPermut3Dpos(:,3),'o','MarkerEdgeColor',[186,212,244]/256);
% plot3(popTarget3Dpos(:,1),popTarget3Dpos(:,2),popTarget3Dpos(:,3),'.','Color',[0.851,0.325,0.098]);
% if samePop == false
%     plot3(popSource3Dpos(:,1),popSource3Dpos(:,2),popSource3Dpos(:,3),'.','Color',[0,0.498,0]);
% end
% axis equal


%% Effectuate the simulations of the random permutations of the popTarget
[dnSimu, ~] = simulateSpatialDisp(NNExp, pops, rowPermut, samePop, nTarget, PARAMS);

%% Calculate exp and simulations cdf and their dispersions individualy
[expCDFs, simuCDFs] = formatCdfs(dnExp, dnSimu, nTarget, PARAMS);

%% Display all the CDFs
figure
displayCDFs(expCDFs, simuCDFs, PARAMS)

% % 
% % figure
% % for simu = 1:PARAMS.numPermut % each simulation is treated individually
% %     subplot(2,3,simu)
% %     hold on
% %     [fSimu(:,simu),xSimu(:,simu),fSimuLo5(:,simu),fSimuUp5(:,simu)] = ecdf(dnSimu(:,simu),'alpha',0.05);
% %     [~,~,fSimuLo1(:,simu),fSimuUp1(:,simu)] = ecdf(dnSimu(:,simu),'alpha',0.01);
% %     plot(xSimu(:,simu),fSimu(:,simu),'.-');
% %     plot(xSimu(:,simu),fSimuLo5(:,simu), '.--');
% %     plot(xSimu(:,simu),fSimuUp5(:,simu), '.--');
% % end
% 
% % 
% % subplot(2,3,6)
% % plot(xSimus,fSimus,'.-');
% % hold on
% % plot(xSimus,[fSimusLo5 fSimusUp5], '.--');
% 
% 8.0
% % Display figure
% figTitle = ['Type II (random permutations in I and II) CellDiameter=',num2str(CellDiameter,2),'{\mu}m'];
% figSavePath = [path,name,'Case',num2str(k),'_Analysis03Fig1'];
% 
% GrandCdf = displaySimuPPQ(Grand, GExp, figTitle, figSavePath, PARAMS);
% 
% save([path,name,'Case',num2str(k),'_Analysis03NN'],'dnExp','GExp','r','Grand','GrandCdf');
% 
fullResults = {};
fullResults.cellDiameter = CellDiameter;
fullResults.dnExp = dnExp;
fullResults.expCDFs = expCDFs;
fullResults.dnSimu = dnSimu;
fullResults.simuCDFs = simuCDFs;
% fullResults.GExp = GExp;
% fullResults.GrandCdf = GrandCdf;

end

function displayProbMap(NNExp, probMap, PARAMS)
% Display probaMap as a colormap applied on the permutable population

all3Dsource = table2array(NNExp(NNExp.cellType == pops.popSource, {'pos3D'}));
all3Dpermut = table2array(NNExp(rowPermut, {'pos3D'}));

colors = colormap(parula(log(max(probMap))/log(PARAMS.effectStrength)+1));
for bioCell = 1:size(all3Dpermut,1)
    colorScat(bioCell,:) = colors(log(probMap(bioCell))/log(PARAMS.effectStrength)+1,:);
end

figure
% scatter3(all3Dpermut(:,1),all3Dpermut(:,2),all3Dpermut(:,3),[],colors(probMap(:),:));
scatter3(all3Dpermut(:,1),all3Dpermut(:,2),all3Dpermut(:,3),[],colorScat,'.');
hold on
scatter3(all3Dsource(:,1),all3Dsource(:,2),all3Dsource(:,3),'.','MarkerFaceColor',[0.098,0.325,0.851]);

end

function dn = findNN(popSource3Dpos, popTarget3Dpos, samePop, pops)
% Recreate the distance table of the 2 cell populations (all neighbours)
% sorted by size and only takes into account the firsts 2 distances

allN = pdist2(popSource3Dpos,popTarget3Dpos,'euclidean','Smallest', 2);
% extract the nearest neighbour distance
if samePop
    % When the target and source are the same then the minimum distance is always 0
    dn = allN(2,:);
else
    % When the target and source are the different
    dn = allN(1,:);
end
if min(dn) == 0
    fprintf('WARNING: nearest neighbour for exp t%d vs t%d is null\n',pops.popSource,pops.popTarget)
end

end

function [dnSimu, Grand] = simulateSpatialDisp(NNExp, pops, rowPermut, samePop, nTarget, PARAMS)
% Simulates random draws of the popTarget population among the popPermut
% population, considering the effect of the popSource population on the
% next draw.

dnSimu = zeros(nTarget,PARAMS.numPermut);

% Calculate all the cell-cell distances to use in the spatial effect
% simulation effect
cell2CellDist = pdist2(table2array(NNExp(:,{'pos3D'})),table2array(NNExp(:,{'pos3D'})));


for perm = 1:PARAMS.numPermut % for each permutation run
    if mod(perm,100) == 0
        fprintf('Running permutation %d\n',perm);
    end
    % Initiate a probability map => prob = 1 for the permutable population and
    % 0 for the rest of the population
    probMap = double(rowPermut);
    popTargetSimu{perm} = table;

    if samePop 
        % Requires adaptation of the popSource at every step + an adaptation at
        % every step of the draw
        for bioCell = 1:nTarget
            tempCell = datasample(NNExp,1,'Weights',probMap);
            % Take cell out of the permutable population
            % be careful at not putting it back at more than 0 by an attractive or repulsive effect !
            probMap(tempCell.cellID) = 0;
            
            % Adapt the probability map of the other cells based on the draw.
            probMap = adaptProbMap(tempCell.cellID, probMap, cell2CellDist, PARAMS);

         
            % Add the drawn cell to the complete simulated draw
            popTargetSimu{perm} = [popTargetSimu{perm};tempCell];
        end
        % In this case target and source populations are the same
        popSourceSimu{perm} = popTargetSimu{perm};
        
    else
        % Adapt the probability map once and for all (no adaptation if there is no
        % spatial effect) => Could also be done once and for all at the beginning 
        % of the simulations
        
        % extract source cells IDs 
        sourceCells = table2array(NNExp(NNExp.cellType == pops.popSource, {'cellID'}));            
        probMap = adaptProbMap(sourceCells, probMap, cell2CellDist, PARAMS);
        
        % Then make a single draw without replacing the cells
        popTargetSimu{perm} = datasample(NNExp,nTarget,'Replace',false,'Weights',probMap);
        popSourceSimu{perm} = NNExp(NNExp.cellType == pops.popSource, :);
    end
    
    % Find the Nearest neighbours in the population
    dnSimu(:,perm) = findNN(popSourceSimu{perm}.pos3D, popTargetSimu{perm}.pos3D, samePop, pops);
    
    % Calculate the cdf histogram
    for i=1:size(PARAMS.binSize,2)
        Grand(perm,i)= sum(dnSimu(:,perm)<PARAMS.binSize(i)) / numel(dnSimu(:,perm));
    end
end

end

function newProbMap = adaptProbMap(sourceCells, oriProbMap, cell2CellDist, PARAMS)
% Adapt the draw probability map based on the desired effect, the drawned
% cell(s)
% SourceCells = cellID of the source population to integrate into the
% probMap
% -> 1 in the case of on the fly adaptation (same source and target pops)
% -> N in the case of static probMap (different source and target pops)

tempProbMap = oriProbMap;

for bioCell = 1:length(sourceCells)
    distsN = cell2CellDist(:,sourceCells(bioCell));
    affectedCells = distsN<PARAMS.effectRange;
    tempProbMap(affectedCells) = tempProbMap(affectedCells)*PARAMS.effectStrength;
end

% if later choose to change the multiplication factor by an additive term,
% multiply tempProMap by a vector containing the 0 positions of the
% oriProbMap (1 otherwise) to ensure conservation of the forbidden
% positions
newProbMap = tempProbMap;

end

function [expCDFs, simuCDFs] = formatCdfs(dnExp, dnSimu, nTarget, PARAMS)

% Calculate the CDFs of the experimental population
[expCDFs.f,expCDFs.x,expCDFs.f5,expCDFs.f95] = ecdf(dnExp,'alpha',0.05);
[~,~,expCDFs.f1,expCDFs.f99] = ecdf(dnExp,'alpha',0.01);

% %% unconclusive version using the greenwood estimator of the variance for each individual cdf
% % Calculate the CDFs of the simulated populations
% for simu = 1:PARAMS.numPermut % each simulation is treated individually
%     % WARNING: Sometimes if popSource = popTarget, the nearest neighbour are
%     % symetrical and return the same distance => are counted as one step in the
%     % cdf of twice the size which messes up the number of positions => error in
%     % table and error in averaging, corrected by interpolation
%     [simuCDFs.indiv{simu}.f,simuCDFs.indiv{simu}.x,simuCDFs.indiv{simu}.fsLo5,simuCDFs.indiv{simu}.fsUp5] = ...
%         ecdf(dnSimu(:,simu),'alpha',0.05);
%     [~,~,simuCDFs.indiv{simu}.fsLo1,simuCDFs.indiv{simu}.fsUp1] = ecdf(dnSimu(:,simu),'alpha',0.01);
%     
%     
%     % In order to avoid dimension mismatch, cdfs are interpolated on fixed
%     % abscissa with fixed periodicity before being merged together
%     simuCDFs.xs(:,simu) = PARAMS.binSize;
%     simuCDFs.fs(:,simu) = interp1([0;unique(simuCDFs.indiv{simu}.x)],...
%         simuCDFs.indiv{simu}.f,PARAMS.binSize);
%     simuCDFs.fsLo5(:,simu) = interp1([0;unique(simuCDFs.indiv{simu}.x)],...
%          simuCDFs.indiv{simu}.fsLo5,PARAMS.binSize);
%     simuCDFs.fsUp5(:,simu) = interp1([0;unique(simuCDFs.indiv{simu}.x)],...
%         simuCDFs.indiv{simu}.fsUp5,PARAMS.binSize);
%     simuCDFs.fsLo1(:,simu) = interp1([0;unique(simuCDFs.indiv{simu}.x)],...
%         simuCDFs.indiv{simu}.fsLo1,PARAMS.binSize);
%     simuCDFs.fsUp1(:,simu) = interp1([0;unique(simuCDFs.indiv{simu}.x)],...
%         simuCDFs.indiv{simu}.fsUp1,PARAMS.binSize);
% end
% 
% % Switch NaNs for 1 (for high values of x, the interpolation returns NaN
% % when there are no available values
% simuCDFs.fs(isnan(simuCDFs.fs)) = 1;
% % simuCDFs.fLo5 = 1 % ENVELOPES SHOULD NOT BE SUBJECT TO THE SAME SWITCH !!!
% % simuCDFs.fUp5 = 1 % ENVELOPES SHOULD NOT BE SUBJECT TO THE SAME SWITCH !!!
% % simuCDFs.fLo1 = 1 % ENVELOPES SHOULD NOT BE SUBJECT TO THE SAME SWITCH !!!
% % simuCDFs.fUp1 = 1 % ENVELOPES SHOULD NOT BE SUBJECT TO THE SAME SWITCH !!!


%% New version using a "simpler" and more manual percentile approach (still kept 
% Greenwood for the experimental cdf. (same process as original WS's)

% Calculate the CDFs of the simulated populations
for simu = 1:PARAMS.numPermut % each simulation is treated individually
    [simuCDFs.indiv{simu}.f,simuCDFs.indiv{simu}.x] = ecdf(dnSimu(:,simu));
    
    % In order to avoid dimension mismatch, cdfs are interpolated on fixed
    % abscissa with fixed periodicity before being merged together
    %     simuCDFs.xs(:,simu) = PARAMS.binSize;
    simuCDFs.fs(:,simu) = interp1([0;unique(simuCDFs.indiv{simu}.x)],...
        [0;simuCDFs.indiv{simu}.f(2:end)],PARAMS.binSize);
end

% Calculate the median simulation
simuCDFs.fs(isnan(simuCDFs.fs)) = 1;
simuCDFs.x = PARAMS.binSize;
% Use percentiles of the individually simulated cdfs to define an envelope
tempOutput = prctile(simuCDFs.fs,[1 5 25 50 75 95 99],2);
simuCDFs.f1pc = tempOutput(:,1);
simuCDFs.f5pc = tempOutput(:,2);
simuCDFs.f25pc = tempOutput(:,3);
simuCDFs.f50pc = tempOutput(:,4);
simuCDFs.f75pc = tempOutput(:,5);
simuCDFs.f95pc = tempOutput(:,6);
simuCDFs.f99pc = tempOutput(:,7);
% Add additionnal measurements
simuCDFs.fmean = mean(simuCDFs.fs,2);
simuCDFs.fstd = std(simuCDFs.fs,[],2);

end

function displayCDFs(expCDFs, simuCDFs, PARAMS)

ylabel('Cumulative cell frequency');
xlabel('Distance to nearest neighbor (µm)');

% Use 2 colors only
colors = [lines(2) [0.5;0.5]];

hold on
% Plot the experimental data
h(1) = plot(expCDFs.x,expCDFs.f,'linewidth',2,'Color',colors(2,1:3));
h(2) = plot(expCDFs.x,expCDFs.f5,'linewidth',1,'Color',colors(2,:));
plot(expCDFs.x,expCDFs.f95,'linewidth',1,'Color',colors(2,:));
h(3) = plot(expCDFs.x,expCDFs.f1,'--','linewidth',1,'Color',colors(2,:));
plot(expCDFs.x,expCDFs.f99,'--','linewidth',1,'Color',colors(2,:));

% Plot the simulation dispersions
h(4) = plot(simuCDFs.x,simuCDFs.f50pc,'linewidth',2,'color',colors(1,1:3));
h(5) = plot(simuCDFs.x,simuCDFs.f5pc,'linewidth',1,'color',colors(1,:));
plot(simuCDFs.x,simuCDFs.f95pc,'linewidth',1,'color',colors(1,:));
h(6) = plot(simuCDFs.x,simuCDFs.f1pc,'--','linewidth',1,'color',colors(1,:));
plot(simuCDFs.x,simuCDFs.f99pc,'--','linewidth',1,'color',colors(1,:));

text(0.5,0.95,'95% and 99% intervals')
text(0.5,0.9,[num2str(PARAMS.numPermut),' random perm.'])

% force the axes
axis(PARAMS.axis);

% display legend
legend(h,{'Experimental data','95% enveloppe','99% enveloppe', ...
    regexprep(sprintf('Simulated data (%s)',PARAMS.model),'_',' '),'95% enveloppe',...
    '99% enveloppe'},'Location','southeast');
legend boxoff;

end





