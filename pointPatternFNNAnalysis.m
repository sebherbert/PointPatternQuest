


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
PARAMS.cellDiameter = mean(NNExp.nearestNeighbour);


%% Extract the 3D positions of the 2 cell populations of interest
popSource3Dpos = table2array(NNExp(NNExp.cellType == pops.popSource, {'pos3D'}));

popTarget3Dpos = table2array(NNExp(NNExp.cellType == pops.popTarget, {'pos3D'}));
nTarget = length(popTarget3Dpos);

rowPermut = logical(sum(NNExp.cellType == pops.popPermut,2));

if pops.popTarget == pops.popSource
    PARAMS.samePop = true;
else 
    PARAMS.samePop = false;
end

% Find nearest neighbour
dnExp = findNN(popSource3Dpos, popTarget3Dpos, PARAMS.samePop, pops)';

% If the range of the model has to be fitted to the cdf50 of the NN
% distribution
if PARAMS.useRangeCDF50 == 1 % Boolean to replace the fit RANGE initial value based on the cdf
    PARAMS.optiR0 = double(median(dnExp));
end

% Interpolate the CDF on a fixed scale and estimate the envelopes
expCDFs = formatCdfsExp(dnExp, PARAMS);

%% Use an optimisation function to find the most adapted strength and range parameters 
if PARAMS.doOptimizePar
    bestParams = optimizeParamsCall(NNExp, pops, rowPermut, nTarget, expCDFs,...
        PARAMS);
    fprintf('bestParams: Range = %0.1fµm ; Strength = %0.2f\n',bestParams(1),bestParams(2));
    
    % Use the optimized values to recalculate the simulation
    effect.Range(1) = bestParams(1);
    effect.Strength = bestParams(2);
    [dnSimu, ~] = simulateSpatialDisp(effect, NNExp, pops, rowPermut, nTarget, PARAMS);
    
    simuCDFs = formatCdfsSimu(dnSimu, PARAMS);   
    % Display simulation + experimental
    figure
    displayCDFs(expCDFs, simuCDFs, PARAMS, PARAMS.useRMSMaxDist);
    
    % recalculate RMS and GOF
    diffCdf = expCDFs.fFix' - simuCDFs.f50pc;
    diffCdf(isnan(diffCdf)) = 0;
    medRMS = median(rms(diffCdf));
    text(60,0.5,sprintf('Fitted parameters:\nRange = %0.1fµm\nStrength = %0.2f\nRMS = %0.04f',...
        bestParams(1),bestParams(2),medRMS));
    saveas(gcf,sprintf('%s_fitted model', fullPath));
    saveas(gcf,sprintf('%s_fitted model.png', fullPath));

else % Use preset values
    %% Effectuate the simulations of the random permutations of the popTarget
    effect.Strength = PARAMS.effectStrength;
    effect.Range = PARAMS.effectRange;
    [dnSimu, ~] = simulateSpatialDisp(effect, NNExp, pops, rowPermut, nTarget, PARAMS);
    
    %% Calculate exp and simulations cdf and their dispersions individualy
    simuCDFs = formatCdfsSimu(dnSimu, PARAMS);
    
    %% Display all the CDFs
    if PARAMS.displayIndivCDF
        figure
        displayCDFs(expCDFs, simuCDFs, PARAMS)
    end
    
    medRMS = evaluateRMS(expCDFs,simuCDFs,PARAMS);
    
end

% Structure output
fullResults = {};
fullResults.dnExp = dnExp;
fullResults.expCDFs = expCDFs;
fullResults.dnSimu = dnSimu;
fullResults.simuCDFs = simuCDFs;
if PARAMS.doOptimizePar
    fullResults.fit.Range = bestParams(1);
    fullResults.fit.Strength = bestParams(2);
    fullResults.fit.medRMS = medRMS;
else
    fullResults.medRMS = medRMS;
end
fullResults.PARAMS = PARAMS;

end

function medRMS = evaluateRMS(expCDFs,simuCDFs,PARAMS)
% Calculate the RMS between experimental and simulated curve f50pc

diffCdf = expCDFs.fFix' - simuCDFs.f50pc;
diffCdf(isnan(diffCdf)) = 0;

if PARAMS.useRMSMaxDist ==1
    % only take into account the distance ditribution up to a maximum distance
    % (as a function of the cell diameter)
    medRMS = rms(diffCdf(simuCDFs.x<PARAMS.maxDistFactor*PARAMS.cellDiameter));
else
    medRMS = rms(diffCdf);
end

end


function [bestParams, finalRMS] = optimizeParamsCall(NNExp, pops, rowPermut, nTarget, expCDFs, PARAMS)
% Launch the optimization function for the fit of the strength of range of
% the distribution effect
x0 = [PARAMS.optiR0,PARAMS.optiS0];

effect.Range = x0(1);
effect.Strength = x0(2);
[dnSimu, ~] = simulateSpatialDisp(effect, NNExp, pops, rowPermut, nTarget, PARAMS);
simuCDFs = formatCdfsSimu(dnSimu, PARAMS);

if PARAMS.doDisplayLiveFit
    % Display the background image (exp value)
    colors = lines(2);
    figure
    h(1) = plot(expCDFs.x,expCDFs.f,'linewidth',2,'Color',colors(2,:));
    hold on
    
    % And original fit (could be added to the )
    h(2) = plot(simuCDFs.x,simuCDFs.f50pc,'linewidth',2,'Color',colors(1,:));
    
    ylabel('Cumulative cell frequency');
    xlabel('Distance to nearest neighbor (µm)');
    
    pause(0.1)
end

% options = optimset('PlotFcns',@optimplotfval);
options = optimset('MaxIter',PARAMS.fitMaxIter);

if PARAMS.doDisplayLiveFit
    [bestParams, finalRMS]= fminsearch(@two_varModel, x0, options, NNExp, pops, rowPermut, nTarget, expCDFs, h, PARAMS);
else
    [bestParams, finalRMS]= fminsearch(@two_varModel, x0, options, NNExp, pops, rowPermut, nTarget, expCDFs, [], PARAMS);
end


end


function [medRMS, h] = two_varModel(x0, NNExp, pops, rowPermut, nTarget, expCDFs, h, PARAMS)
% provide the ks test comparing the experimental and simulated
% distributions

effect.Range = x0(1);
effect.Strength = x0(2);

if (effect.Range <= PARAMS.minFitRange || effect.Strength <= PARAMS.minFitStrength)
    medRMS = 100;
    return
end

[dnSimu, ~] = simulateSpatialDisp(effect, NNExp, pops, rowPermut, nTarget, PARAMS);

simuCDFs = formatCdfsSimu(dnSimu, PARAMS);

medRMS = evaluateRMS(expCDFs,simuCDFs,PARAMS);

if PARAMS.doDisplayLiveFit
    % update image
    set(h(2),'YData',simuCDFs.f50pc,'Color',lines(1));
    pause(0.1)
end

fprintf('bestParams: Range = %0.1fµm ; Strength = %0.2f ; RMS = %0.5f\n',effect.Range,effect.Strength,medRMS);

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

function [dnSimu, Grand] = simulateSpatialDisp(effect, NNExp, pops, rowPermut, nTarget, PARAMS)
% Simulates random draws of the popTarget population among the popPermut
% population, considering the effect of the popSource population on the
% next draw.

dnSimu = zeros(nTarget,PARAMS.numPermut);

% Calculate all the cell-cell distances to use in the spatial effect
% simulation effect
cell2CellDist = pdist2(table2array(NNExp(:,{'pos3D'})),table2array(NNExp(:,{'pos3D'})));

Grand = zeros(PARAMS.numPermut, size(PARAMS.binSize,2));

parfor perm = 1:PARAMS.numPermut % for each permutation run
    %     if mod(perm,1000) == 0
    %         fprintf('Running permutation %d\n',perm);
    %     end
    % Initiate a probability map => prob = 1 for the permutable population and
    % 0 for the rest of the population
    probMap = double(rowPermut);
    popTargetSimu{perm} = table;
    
    if PARAMS.samePop
        % Requires adaptation of the popSource at every step + an adaptation at
        % every step of the draw
        for bioCell = 1:nTarget
            tempCell = datasample(NNExp,1,'Weights',probMap);
            % Take cell out of the permutable population
            % be careful at not putting it back at more than 0 by an attractive or repulsive effect !
            probMap(tempCell.cellID) = 0;
            
            % Adapt the probability map of the other cells based on the draw.
            probMap = adaptProbMap(tempCell.cellID, probMap, cell2CellDist, effect);

         
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
        probMap = adaptProbMap(sourceCells, probMap, cell2CellDist, effect);
        
        % Then make a single draw without replacing the cells
        popTargetSimu{perm} = datasample(NNExp,nTarget,'Replace',false,'Weights',probMap);
        popSourceSimu{perm} = NNExp(NNExp.cellType == pops.popSource, :);
    end
    
    % Find the Nearest neighbours in the population
    dnSimu(:,perm) = findNN(popSourceSimu{perm}.pos3D, popTargetSimu{perm}.pos3D, PARAMS.samePop, pops);
    
    %     % Calculate the cdf histogram
    %     for i=1:size(PARAMS.binSize,2)
    %         Grand(perm,i)= sum(dnSimu(:,perm)<PARAMS.binSize(i)) / numel(dnSimu(:,perm));
    %     end
end

end

function newProbMap = adaptProbMap(sourceCells, oriProbMap, cell2CellDist, effect)
% Adapt the draw probability map based on the desired effect, the drawned
% cell(s)
% SourceCells = cellID of the source population to integrate into the
% probMap
% -> 1 in the case of on the fly adaptation (same source and target pops)
% -> N in the case of static probMap (different source and target pops)

tempProbMap = oriProbMap;

for bioCell = 1:length(sourceCells)
    distsN = cell2CellDist(:,sourceCells(bioCell));
    affectedCells = distsN<effect.Range;
    tempProbMap(affectedCells) = tempProbMap(affectedCells)*effect.Strength;
end

% if later choose to change the multiplication factor by an additive term,
% multiply tempProMap by a vector containing the 0 positions of the
% oriProbMap (1 otherwise) to ensure conservation of the forbidden
% positions
newProbMap = tempProbMap;

end

function expCDFs = formatCdfsExp(dnExp, PARAMS)
% Interpolate the CDF on a fixed scale and estimate the envelopes

% Calculate the CDFs of the experimental population
[expCDFs.f,expCDFs.x,expCDFs.f5,expCDFs.f95] = ecdf(dnExp,'alpha',0.05);
[~,~,expCDFs.f1,expCDFs.f99] = ecdf(dnExp,'alpha',0.01);

% Interpolate the CDF on a fixed scale
expCDFs.fFix = interp1([0;unique(expCDFs.x)],...
    [0;expCDFs.f(2:end)],PARAMS.binSize);
expCDFs.xFix = PARAMS.binSize;
end

function simuCDFs = formatCdfsSimu(dnSimu, PARAMS)

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

function displayCDFs(expCDFs, simuCDFs, PARAMS, dispMaxCdf)

ylabel('Cumulative cell frequency');
xlabel('Distance to nearest neighbor (µm)');

% Use 2 colors only (3rd for the max cdf)
colors = [lines(3) [0.5;0.5;1]];

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

% text(0.5,0.95,'95% and 99% intervals')
% text(0.5,0.9,[num2str(PARAMS.numPermut),' random perm.'])

% Set up the legend
legs = {'Experimental data','95% enveloppe','99% enveloppe', ...
    regexprep(sprintf('Simulated data (%s)',PARAMS.model),'_',' '),'95% enveloppe',...
    '99% enveloppe'};

if dispMaxCdf
    % last value used for RMS calculation
    xPos = (sum(simuCDFs.x<PARAMS.maxDistFactor*PARAMS.cellDiameter));
    h(7) = plot(simuCDFs.x(xPos),simuCDFs.f50pc(xPos),'o',...
        'MarkerSize',10, 'MarkerEdgeColor',[0,0.5,0],...
        'MarkerFaceColor',[0,0.5,0]);
    % Adapt the legend
    legs = [legs sprintf('maxCDF = %0.1fum/%dcellDia',PARAMS.maxDistFactor*PARAMS.cellDiameter,...
        PARAMS.maxDistFactor)];
end

% force the axes
axis(PARAMS.axis);

% display legend
legend(h,legs,'Location','southeast');
legend boxoff;

end





