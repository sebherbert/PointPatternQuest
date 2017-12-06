

function combineAndCompareNN()
% Sets up automated ks tests in between the different populations
% Samples are fixed at 3 populations (young, old, drug) with 3 samples in
% each population.
% For each individual there are 3 tests (II on II, III on II, III on III)
% and 4 possible regions on the brain => all, da, dl and dm defined by the
% variable k = 1,2,3,4 respectively. In pratical case, the variable 4 is
% likely not going to be used.

%% Load up the samples per populations
% The script can handle a different number of individuals and incomplete
% conditions per individual.

close all

PARAMS = {};
PARAMS.axes = [0 10 0 1];

%% Initialize the parameters
% inds = {'ind1','ind2','ind3'}; % names of the individual fishes => Number
% of individual per condition can change => Not required anymore
bps = {'allSample','dm','da'}; % parts of the brain (not da for the moment since not tested in a brain)
NNtests = {'t2vst2','t3vst3','t3vst2'}; % tested cell types
conds = {'Young','Old','Drug'}; % Conditions to be tested => CHECK THE ORDER WITH allDATA STRUCTURE

% temp = load('/media/sherbert/Data/Projects/OG_projects/Project6_ND_distribPattern/dataFolder/loader.mat'); % temporary

% tempFold = uipickfiles('Prompt','Please, select the output folder');
tempFold = {'/media/sherbert/Data/Projects/OG_projects/Project6_ND_distribPattern/dataFolder/test'};
PARAMS.outputFold = [tempFold{1} filesep];


% Load populations
allData = {};
for condition = 1:numel(conds) % for each condition
    multiFoldPath = uipickfiles('Prompt',... % temporarilly killed
        sprintf('Select the folders of individual in the %s population',conds{condition}));
    %     multiFoldPath = temp.foldPaths{condition};
    for folderOI = 1:length(multiFoldPath) % for each selected folder
        folderInfo = dir(multiFoldPath{folderOI}); % find objects into the folder
        for obj = 1:length(folderInfo)
            if ~folderInfo(obj).isdir % if object is not a folder
                for bp = 1:numel(bps)
                    if contains(folderInfo(obj).name,sprintf('_%s',bps{bp}))
                        fprintf('Found file: %s\n',folderInfo(obj).name);
                        dataTemp = load([folderInfo(obj).folder filesep folderInfo(obj).name]);
                        tempField = fieldnames(dataTemp);
                        if numel(tempField) ~= 1
                            fprintf('Found file has an abnormal number of fields');
                        else
                            allData{condition}.(sprintf('ind%d',folderOI)).(bps{bp}) = dataTemp.(tempField{1});
                        end
                        clear dataTemp
                    end
                end
            end
        end
        [~,indName,~] = fileparts(multiFoldPath{folderOI});
        individualNames.(conds{condition}).(sprintf('ind%d',folderOI)) = indName;
        clear indName;
    end
end

%% Save connection of individual names and ind$N
save([PARAMS.outputFold 'individualNames'],'individualNames');


%% Display inter and intra (1 color each) subplot 3(parts)x3tests
lineColors = lines(2); % compare populations 2x2
for condition1 = 1:numel(conds)
    for condition2 = condition1+1:numel(conds)
        conditions = {conds{condition1}, conds{condition2}}; % Conditions to be tested
        displaySubPlots({allData{condition1},allData{condition2}}, conditions, bps, NNtests, lineColors, PARAMS);
    end
end

% Pool all the conditions
lineColors = lines(numel(conds));
displaySubPlots({allData{1},allData{2},allData{3}},conds,bps,NNtests,lineColors,PARAMS);
% 
% for condition = 1:numel(conds) 
%     displayPerCondition(PARAMS)
% end
% 


%% ks tests
% ks intra population
% For each individual condition to be tested
ksIntra = ksTestIntra({allData{1},allData{2},allData{3}},conds,bps,NNtests);

% ks inter population => Automatize
ksInter = {}; % initialize structure
for condition1 = 1:numel(conds)
    for condition2 = condition1+1:numel(conds)
        conditions = {conds{condition1}, conds{condition2}}; % Conditions to be tested
        ksInter = ksTestInter({allData{condition1},allData{condition2}},conditions,bps,NNtests,ksInter);
    end
end

% Reshape for reading in 3x3 tables cf excel sheets (name rows and lines)
rsKSpInter = reshapeKS(ksInter,'ksInter',individualNames);
rsKSpIntra = reshapeKS(ksIntra,'ksIntra',individualNames);

save([PARAMS.outputFold 'ksResults'],'ksInter','ksIntra','rsKSpInter','rsKSpIntra');

end

function rsKSp = reshapeKS(ksStruct,structName,individualNames)
% Reshape the KS tests to make it more readable

rsKSp = {};
crossTests = fieldnames(ksStruct);
for crossTest = 1:numel(crossTests)
    crossInds = fieldnames(ksStruct.(crossTests{crossTest}));
    for crossInd = 1:numel(crossInds)
        localBps = fieldnames(ksStruct.(crossTests{crossTest}).(crossInds{crossInd}));
        for localBp = 1:numel(localBps)
            localNNts = fieldnames(ksStruct.(crossTests{crossTest}).(crossInds{crossInd}). ...
                (localBps{localBp}));
            for localNNt = 1:numel(localNNts)
                
                tempInd = regexp((crossInds{crossInd}),'(?<=ind[^0-9]*)\d*','match');
                if numel(tempInd) > 2
                    fprintf('WARNING: Reshaping of %s failed (too many ind)',structName);
                    break
                elseif numel(tempInd) < 2
                    fprintf('WARNING: Reshaping of %s failed (too few ind)',structName);
                end
                
                % Save KSs in an array
                ksTemp = ksStruct.(crossTests{crossTest}).(crossInds{crossInd}). ...
                    (localBps{localBp}).(localNNts{localNNt}).p;
                
                rsKSp.(crossTests{crossTest}).(localNNts{localNNt}). ...
                    (localBps{localBp})(str2double(tempInd{1}),str2double(tempInd{2})) = ...
                    ksTemp;

                
                % must keep track of every individual in the table to recover the proper
                % names afterward => here we save the ind#
                if isfield(rsKSp.(crossTests{crossTest}).(localNNts{localNNt}), ([localBps{localBp} '_AllRowNames']))
                    rsKSp.(crossTests{crossTest}).(localNNts{localNNt}).([localBps{localBp} '_AllRowNames'])(end+1) = ...
                        str2double(tempInd{1});
                    %                         horzcat(rsKSp.(crossTests{crossTest}).(localNNts{localNNt}).([localBps{localBp} '_AllRowNames']),...
                    rsKSp.(crossTests{crossTest}).(localNNts{localNNt}).([localBps{localBp} '_AllVarNames'])(end+1) = ...
                        str2double(tempInd{2});
                    %                         horzcat(rsKSp.(crossTests{crossTest}).(localNNts{localNNt}).([localBps{localBp} '_AllVarNames']),...

                else
                    rsKSp.(crossTests{crossTest}).(localNNts{localNNt}).([localBps{localBp} '_AllRowNames']) = ...
                        str2double(tempInd{1});
                    rsKSp.(crossTests{crossTest}).(localNNts{localNNt}).([localBps{localBp} '_AllVarNames']) = ...
                        str2double(tempInd{2});
                end
                    
            end
        end
    end
end

% Change arrays for tables with individual names
crossTests = fieldnames(rsKSp);
for crossTest = 1:numel(crossTests)
    localNNts = fieldnames(rsKSp.(crossTests{crossTest}));
    for localNNt = 1:numel(localNNts)
        localBps = fieldnames(rsKSp.(crossTests{crossTest}).(localNNts{localNNt}));
        for localBp = 1:numel(localBps)
            
            % check if the field name is filled with data or legends
            if contains(localBps{localBp},'Names')
                continue
            end
            
            tempTable = array2table(rsKSp.(crossTests{crossTest}).(localNNts{localNNt}).(localBps{localBp}));
            % Create the names array
            localConds = regexp(crossTests{crossTest},'vs','split');
                        
            % Set the row names
            clear allRN
            UniqRowNames = unique(rsKSp.(crossTests{crossTest}).(localNNts{localNNt}).([localBps{localBp} '_AllRowNames']));
            for indN = 1:size(tempTable,1)
                allRN{indN} = individualNames.(localConds{1}).(sprintf('ind%d', UniqRowNames(indN)));
                if ~isvarname(allRN{indN})
                    allRN{indN} = sprintf('Ind_%s',allRN{indN});
                end
            end
            tempTable.Properties.RowNames = allRN;

            % Set the column names
            clear allVN
            if numel(localConds)==1
                allVN = allRN;
            else
                UniqVarNames = unique(rsKSp.(crossTests{crossTest}).(localNNts{localNNt}).([localBps{localBp} '_AllVarNames']));
                for indN = 1:size(tempTable,2)
                    allVN{indN} = individualNames.(localConds{2}).(sprintf('ind%d', UniqVarNames(indN)));
                    if ~isvarname(allVN{indN})
                        allVN{indN} = sprintf('Ind_%s',allVN{indN});
                    end
                end
            end
            tempTable.Properties.VariableNames = allVN;
            
 
            rsKSp.(crossTests{crossTest}).(localNNts{localNNt}).(localBps{localBp}) = tempTable;
        end
    end
end

end

function ksIntra = ksTestIntra(fullData,conditions,bps,NNtests)
% Calculates and return the individual and average ks tests inside a single
% population

% ks intra pop
ksIntra = {};
% example => ksIntra.condition.ind1.bps.NNtests = kstest2;
for NNtest = 1:numel(NNtests) % which test
    for bp = 1:numel(bps) % which brain part
        for condition = 1:numel(conditions) % which condition
            % list every individual in condition
            indNames = fieldnames(fullData{condition});
            for ind1 = 1:numel(indNames) % which individual 1 
                ind1Name = indNames{ind1};
                for ind2 = 1:numel(indNames) % which individual 2
                    
                    ind2Name = indNames{ind2};
                    
                    % if this brain part is not treated in any of these individual then
                    % skip it
                    if ~isfield( fullData{condition}.(ind1Name) , (bps{bp}))
                        continue
                    elseif ~isfield( fullData{condition}.(ind2Name) , (bps{bp}))
                        continue
                    end
                    
                    dn1 = fullData{condition}.(ind1Name).(bps{bp}).(NNtests{NNtest}).dn;
                    dn2 = fullData{condition}.(ind2Name).(bps{bp}).(NNtests{NNtest}).dn;
                    indsName = sprintf('%svs%s',ind1Name,ind2Name);
                    
                    [h,p,k] = kstest2(dn1,dn2);
                    
                    ksIntra.(conditions{condition}).(indsName).(bps{bp}).(NNtests{NNtest}).h = h;
                    ksIntra.(conditions{condition}).(indsName).(bps{bp}).(NNtests{NNtest}).p = p;
                    ksIntra.(conditions{condition}).(indsName).(bps{bp}).(NNtests{NNtest}).k = k;
                end
            end
        end
    end
end

end

function ksInter = ksTestInter(fullData,conditions,bps,NNtests,ksInter)
% Calculates and return the individual and average ks tests inside a single
% population

% ks inter pop
% example => ksInter.condition.ind1vs1ind2.bps.NNtests = kstest2;
for NNtest = 1:numel(NNtests) % which test
    for bp = 1:numel(bps) % which brain part
        for condition1 = 1:numel(conditions) % which condition1

            % list every individual in condition 1
            indNamesCond1 = fieldnames(fullData{condition1});

            for condition2 = condition1+1:numel(conditions) % which condition2
            
                % list every individual in condition 2
                indNamesCond2 = fieldnames(fullData{condition2});
                
                for ind1 = 1:numel(fieldnames(fullData{condition1})) % which individual1
                    ind1Name = indNamesCond1{ind1};
                    for ind2 = 1:numel(fieldnames(fullData{condition2})) % which individual2
                        ind2Name = indNamesCond2{ind2};
                        
                        % if this brain part is not treated in any of these individual then
                        % skip it
                        if ~isfield( fullData{condition1}.(ind1Name) , (bps{bp}))
                            continue
                        elseif ~isfield( fullData{condition2}.(ind2Name) , (bps{bp}))
                            continue
                        end
                        
                        % Simplify names of the 2 pops of interest
                        dn1 = fullData{condition1}.(ind1Name).(bps{bp}).(NNtests{NNtest}).dn;
                        dn2 = fullData{condition2}.(ind2Name).(bps{bp}).(NNtests{NNtest}).dn;
                        
                        % ks test per se
                        [h,p,k] = kstest2(dn1,dn2);
                        
                        % Allocate in structure
                        condsName = sprintf('%svs%s',conditions{condition1},conditions{condition2});
                        indsName = sprintf('%svs%s',ind1Name,ind2Name);
                        ksInter.(condsName).(indsName).(bps{bp}).(NNtests{NNtest}).h = h;
                        ksInter.(condsName).(indsName).(bps{bp}).(NNtests{NNtest}).p = p;
                        ksInter.(condsName).(indsName).(bps{bp}).(NNtests{NNtest}).k = k;
                    end
                end
            end
        end
    end
end

end

function displaySubPlots(fullData,conditions,bps,NNtests,lineColors,PARAMS)
% Display inter and intra (1 color each) subplot 3(parts of brain)x3 tests
figTitle = sprintf('Conditions: %s', conditions{1});
figSaveName = sprintf('%s', conditions{1});

for nbrConds = 2:numel(conditions)
    figTitle = sprintf('%s vs %s', figTitle, conditions{nbrConds});
    figSaveName = sprintf('%svs%s', figSaveName, conditions{nbrConds});
end
figure('Name',figTitle);

r=0:0.1:14; % bin size for the ecdf for as long as I don't auto

for NNtest = 1:numel(NNtests) % which test
    for bp = 1:numel(bps) % which brain part
        subplot(3,3,(bp-1)*3+NNtest)
        hold on
        hExp = {};
        hSimu = {};
        clegs = [];
        legs = [];
        title(sprintf('NN %s in %s',NNtests{NNtest},bps{bp}));
        for condition = 1:numel(conditions) % which condition
            % list every individual
            indNames = fieldnames(fullData{condition});
            Nind = 0; % total pop in the current test
            for ind = 1:numel(indNames) % which individual
                indName = indNames{ind};
                
                % if this brain part is not treated in this particular individual then
                % skip it
                if ~isfield( fullData{condition}.(indName) , (bps{bp}))
                    continue
                end
                
                % If this statistical test is not treated in this particular individual
                % then skip it
                if ~isfield( fullData{condition}.(indName).(bps{bp}), NNtests{NNtest} )
                    continue
                end
                
                % Number of individual per condition / brainpart / NNtest
                Nind = Nind+1;
                
                % plot the experimental cdf
                dnExp = fullData{condition}.(indName).(bps{bp}).(NNtests{NNtest}).dn;
                %                 dn(1,:) = dn(1,:)+numel(hExp); % for tests only
                [fExp,xExp] = ecdf(dnExp); % calculate associated experimental cdf
                hExp{numel(hExp)+1} = plot(xExp,fExp,'.-','Color',lineColors(condition,:));
                % splot of the theoretical cdf => Should use dnSimu in a later version
                dnSimu = fullData{condition}.(indName).(bps{bp}).(NNtests{NNtest}).GrandCdf.mean;
                %                 dnSimu(1,:) = dnSimu(1,:)+numel(hSimu); % for tests only
                %                 [fSimu,xSimu] = ecdf(dnSimu); % calculate associated
                %                 experimental cdf => for when dnSimu = real
                fSimu = dnSimu;
                xSimu = r;
                hSimu{numel(hSimu)+1} = plot(xSimu,fSimu,'--','Color',[lineColors(condition,:) 0.5]);
            end
            clegs = [clegs, hExp{numel(hExp)}, hSimu{numel(hSimu)}];
            legs = [legs, sprintf('%s (N=%d)',conditions{condition},Nind),...
                strcat(conditions(condition),' simus')];
        end
        legend(clegs,legs,'Location','southeast');
        legend boxoff;
        axis(PARAMS.axes);
        if (((bp-1)*3+NNtest) >= 7)
            xlabel('CDF of NN (cell diam)');
        end
    end
end
      
saveas(gcf,[PARAMS.outputFold figSaveName]);
% saveas(gcf,[PARAMS.outputFold figSaveName],'png'); => Need Full Screen
        
% [f,x] = ecdf(dn);
% hold on
% plot(x,f,'.-');
% plot(r,G,'.-');
% % plot(r,GrandAll.mean);

end
