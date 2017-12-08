


% S=S(S==1|S==2)'; for Analysis 3

function full_Results = pointPatternFNNAnalysis(fullPath, NNExp, pops, PARAMS)
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

% For simpler reading
cellType = NNExp.cellType;
r = PARAMS.binSize;

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

% Normalize distance by cell Diameter
dnExp = dnExp/CellDiameter;

% % calculate the cdf by hand... consider using ecdf for proper method
% for i = 1:size(r,2)
%     GExp(i) = size(dnExp(dnExp<=r(i)),2)/size(dnExp,2);
% end

% Display experimental populations cells => temporary
figure
hold on
plot3(popPermut3Dpos(:,1),popPermut3Dpos(:,2),popPermut3Dpos(:,3),'o');
plot3(popTarget3Dpos(:,1),popTarget3Dpos(:,2),popTarget3Dpos(:,3),'.');
axis equal


%% Effectuate the simulations of the random permutations of the popTarget
[dnSimu, ~] = randomPerm(NNExp, pops, rowPermut, samePop, CellDiameter, nTarget, PARAMS);

%% Calculate exp and simulations cdf and their dispersions individualy
[expCDFs, simuCDFs] = formatCdfs(dnExp, dnSimu, nTarget, PARAMS);

%% Display all the CDFs
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
% 
% % Display figure
% figTitle = ['Type II (random permutations in I and II) CellDiameter=',num2str(CellDiameter,2),'{\mu}m'];
% figSavePath = [path,name,'Case',num2str(k),'_Analysis03Fig1'];
% 
% GrandCdf = displaySimuPPQ(Grand, GExp, figTitle, figSavePath, PARAMS);
% 
% save([path,name,'Case',num2str(k),'_Analysis03NN'],'dnExp','GExp','r','Grand','GrandCdf');
% 
% fullResults = {};
% fullResults.dnExp = dnExp;
% fullResults.GExp = GExp;
% fullResults.GrandCdf = GrandCdf;

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

function [dnSimu, Grand] = randomPerm(NNExp, pops, rowPermut, samePop, CellDiameter, nTarget, PARAMS)
% Simulates random draws of the popTarget population among the popPermut
% population, considering the effect of the popSource population on the
% next draw.

dnSimu = zeros(nTarget,PARAMS.numPermut);

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
            if strcmp(PARAMS.effect,'Repulsion')
                fprintf('ERROR: Spatial effect not ready%s\n', PARAMS.effect);
                probMap = sum(NNExp.cellType == pops.popPermut,2); % TBChanged !
            elseif strcmp(PARAMS.effect,'Attraction')
                fprintf('ERROR: Spatial effect not ready%s\n', PARAMS.effect);
                probMap = sum(NNExp.cellType == pops.popPermut,2); % TBChanged !
            elseif ~strcmp(PARAMS.effect,'None')
                fprintf('ERROR: Unknown spatial effect %s\n', PARAMS.effect);
                break
            end
            
            % Add the drawn cell to the complete simulated draw
            popTargetSimu{perm} = [popTargetSimu{perm};tempCell];
        end
        
        % In this case target and source populations are the same
        popSourceSimu{perm} = popTargetSimu{perm};
        
    else
        % Adapt the probability map once and for all (no adaptation if there is no
        % spatial effect)        
        if strcmp(PARAMS.effect,'Repulsion')
            fprintf('ERROR: Spatial effect not ready%s\n', PARAMS.effect);
            probMap = sum(NNExp.cellType == pops.popPermut,2); % TBChanged ! 
        elseif strcmp(PARAMS.effect,'Attraction')
            fprintf('ERROR: Spatial effect not ready%s\n', PARAMS.effect);
            probMap = sum(NNExp.cellType == pops.popPermut,2); % TBChanged ! 
        elseif ~strcmp(PARAMS.effect,'None')
            fprintf('ERROR: Unknown spatial effect %s\n', PARAMS.effect);
            break
        end
        
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

% Normalize distance by cell Diameter
% Grand = Grand/CellDiameter;
dnSimu = dnSimu/CellDiameter;

end

function [expCDFs, simuCDFs] = formatCdfs(dnExp, dnSimu, nTarget, PARAMS)

% Calculate the CDFs of the experimental population
[expCDFs.f,expCDFs.x,expCDFs.fLo5,expCDFs.fUp5] = ecdf(dnExp,'alpha',0.05);
[~,~,expCDFs.fLo1,expCDFs.fUp1] = ecdf(dnExp,'alpha',0.01);

% Calculate the CDFs of the simulated populations

% clear fSimuLo1 fSimuUp1 fSimuLo5 fSimuUp5 xSimu fSimu
simuCDFs.fs = zeros(nTarget+1,PARAMS.numPermut);
simuCDFs.xs = zeros(nTarget+1,PARAMS.numPermut);
simuCDFs.fsLo5 = zeros(nTarget+1,PARAMS.numPermut);
simuCDFs.fsUp5 = zeros(nTarget+1,PARAMS.numPermut);
simuCDFs.fsLo1 = zeros(nTarget+1,PARAMS.numPermut);
simuCDFs.fsUp1 = zeros(nTarget+1,PARAMS.numPermut);
for simu = 1:PARAMS.numPermut % each simulation is treated individually
    % ERROR: Sometimes if popSource = popTarget, the nearest neighbour are
    % symetrical and return the same distance => are counted as one step in the
    % cdf of twice the size which messes up the number of positions => error in
    % table and error in averaging
    [simuCDFs.fs(:,simu),simuCDFs.xs(:,simu),simuCDFs.fsLo5(:,simu),simuCDFs.fsUp5(:,simu)] = ecdf(dnSimu(:,simu),'alpha',0.05);
    [~,~,simuCDFs.fsLo1(:,simu),simuCDFs.fsUp1(:,simu)] = ecdf(dnSimu(:,simu),'alpha',0.01);
end

% Calculate the average simulation
simuCDFs.f = mean(simuCDFs.fs,2);
simuCDFs.x = mean(simuCDFs.xs,2);
simuCDFs.fLo5 = mean(simuCDFs.fsLo5,2);
simuCDFs.fUp5 = mean(simuCDFs.fsUp5,2);
simuCDFs.fLo1 = mean(simuCDFs.fsLo1,2);
simuCDFs.fUp1 = mean(simuCDFs.fsUp1,2);

end

function displayCDFs(expCDFs, simuCDFs, PARAMS)

figure

ylabel('Cumulative cell frequency');
xlabel('Distance to nearest neighbor (in cell diameter)');

% Use 2 colors only
colors = [lines(2) [0.5;0.5]];

hold on
% Plot the experimental data
h(1) = plot(expCDFs.x,expCDFs.f,'linewidth',2,'Color',colors(2,1:3));
h(2) = plot(expCDFs.x,expCDFs.fLo5,'linewidth',1,'Color',colors(2,:));
plot(expCDFs.x,expCDFs.fUp5,'linewidth',1,'Color',colors(2,:));
h(3) = plot(expCDFs.x,expCDFs.fLo1,'--','linewidth',1,'Color',colors(2,:));
plot(expCDFs.x,expCDFs.fUp1,'--','linewidth',1,'Color',colors(2,:));

% Plot the simulation dispersions

h(4) = plot(simuCDFs.x,simuCDFs.f,'linewidth',2,'color',colors(1,1:3));
h(5) = plot(simuCDFs.x,simuCDFs.fLo5,'linewidth',1,'color',colors(1,:));
plot(simuCDFs.x,simuCDFs.fUp5,'linewidth',1,'color',colors(1,:));
h(6) = plot(simuCDFs.x,simuCDFs.fLo1,'--','linewidth',1,'color',colors(1,:));
plot(simuCDFs.x,simuCDFs.fUp1,'--','linewidth',1,'color',colors(1,:));

text(0.5,0.95,'95% and 99% intervals')
text(0.5,0.9,[num2str(PARAMS.numPermut),' random perm.'])

% force the axes
axis([0 5 0 1]);

% display legend
legend(h,{'Experimental data','95% enveloppe','99% enveloppe', ...
    'Simulated data','95% enveloppe','99% enveloppe'},'Location','southeast');
legend boxoff;


% if 2*ceil(xScale(min(find(NN>0.9))))<14
%     axis([0 2*ceil(xScale(min(find(NN>0.9)))) 0 1]);
% else
%     axis([0 14 0 1]);
% end

end





