

function combineAndCompareNN()
% Sets up automated ks tests in between the different populations
% Samples are fixed at 3 populations (young, old, drug) with 3 samples in
% each population.
% For each individual there are 3 tests (II on II, III on II, III on III)
% and 4 possible regions on the brain => all, da, dl and dm defined by the
% variable k = 1,2,3,4 respectively. In pratical case, the variable 4 is
% likely not going to be used.

%% Load up the samples per populations
% 3 independent loads or automated parsing?
% => Objective is to create 3 structures with each N samples, then 3 brain
% parts and then 3 tests

close all

PARAMS = {};

%% Initialize the parameters
inds = {'ind1','ind2','ind3'}; % names of the individual fishes
bps = {'allSample','dm'}; % parts of the brain (not da for the moment since not tested in a brain)
NNtests = {'t2vst2','t3vst3','t3vst2'}; % tested cell types
conds = {'Young','Drug','Old'}; % Conditions to be tested => CHECK THE ORDER WITH allDATA STRUCTURE

temp = load('/media/sherbert/Data/Projects/OG_projects/Project6_ND_distribPattern/dataFolder/loader.mat'); % temporary

% tempFold = uipickfiles('Prompt','Please, select the output folder');
tempFold = {'/media/sherbert/Data/Projects/OG_projects/Project6_ND_distribPattern/dataFolder/Output'};
PARAMS.outputFold = [tempFold{1} filesep];


% Load populations
allData = {};
for condition = 1:numel(conds) % for each condition
    %     multiFoldPath = uipickfiles('Prompt',sprintf('Select the folders of
    %     individual in the %s population',conds{condition})); % temporarilly
    %     killed
    multiFoldPath = temp.foldPaths{condition};
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
    end
end


% % => Old test version
% dataYoung.ind1.all = dataCombined;
% dataYoung.ind2.all = dataCombined;
% dataYoung.ind3.all = dataCombined;
% dataOld.ind1.all = dataCombined;
% dataOld.ind2.all = dataCombined;
% dataOld.ind3.all = dataCombined;
% dataDrug.ind1.all = dataCombined;
% dataDrug.ind2.all = dataCombined;
% dataDrug.ind3.all = dataCombined;
% 
% dataYoung.ind1.da = dataCombined;
% dataYoung.ind2.da = dataCombined;
% dataYoung.ind3.da = dataCombined;
% dataOld.ind1.da = dataCombined;
% dataOld.ind2.da = dataCombined;
% dataOld.ind3.da = dataCombined;
% dataDrug.ind1.da = dataCombined;
% dataDrug.ind2.da = dataCombined;
% dataDrug.ind3.da = dataCombined;
% 
% dataYoung.ind1.dl = dataCombined;
% dataYoung.ind2.dl = dataCombined;
% dataYoung.ind3.dl = dataCombined;
% dataOld.ind1.dl = dataCombined;
% dataOld.ind2.dl = dataCombined;
% dataOld.ind3.dl = dataCombined;
% dataDrug.ind1.dl = dataCombined;
% dataDrug.ind2.dl = dataCombined;
% dataDrug.ind3.dl = dataCombined;
% 
% allData = {dataYoung, dataDrug, dataOld};


%% Display inter and intra (1 color each) subplot 3(parts)x3tests
lineColors = lines(2); % compare populations 2x2
for condition1 = 1:numel(conds)
    for condition2 = condition1+1:numel(conds)
        conditions = {conds{condition1}, conds{condition2}}; % Conditions to be tested
        displaySubPlots({allData{condition1},allData{condition2}}, conditions, inds, bps, NNtests, lineColors, PARAMS);
    end
end

% conditions = {'Young','Old'}; % Conditions to be tested
% lineColors = lines(numel(conditions));
% displaySubPlots({dataYoung,dataOld},conditions,inds,bps,NNtests,lineColors); 
% 
% conditions = {'Young','Drug'}; % Conditions to be tested
% lineColors = lines(numel(conditions));
% displaySubPlots({dataYoung,dataDrug},conditions,inds,bps,NNtests,lineColors);
% 
% conditions = {'Drug','Old'}; % Conditions to be tested
% lineColors = lines(numel(conditions));
% displaySubPlots({dataDrug,dataOld},conditions,inds,bps,NNtests,lineColors);

% Pool all the conditions
lineColors = lines(numel(conds));
displaySubPlots({allData{1},allData{2},allData{3}},conds,inds,bps,NNtests,lineColors,PARAMS);


%% ks tests
% ks intra population
% For each individual condition to be tested
ksIntra = ksTestIntra({allData{1},allData{2},allData{3}},conds,inds,bps,NNtests);

% ks inter population => Automatize
ksInter = {}; % initialize structure
for condition1 = 1:numel(conds)
    for condition2 = condition1+1:numel(conds)
        conditions = {conds{condition1}, conds{condition2}}; % Conditions to be tested
        %         displaySubPlots({allData{condition1},allData{condition2}}, conditions,inds,bps,NNtests,lineColors);
        ksInter = ksTestInter({allData{condition1},allData{condition2}},conditions,inds,bps,NNtests,ksInter);
    end
end

% Reshape for reading in 3x3 tables cf excel sheets (name rows and lines)

save([PARAMS.outputFold 'ksResults'],'ksInter','ksIntra');

end

function ksIntra = ksTestIntra(fullData,conditions,inds,bps,NNtests)
% Calculates and return the individual and average ks tests inside a single
% population

% ks intra pop
ksIntra = {};
% example => ksIntra.condition.ind1.bps.NNtests = kstest2;
for NNtest = 1:numel(NNtests) % which test
    for bp = 1:numel(bps) % which brain part
        for condition = 1:numel(conditions) % which condition
            for ind1 = 1:numel(fieldnames(fullData{condition})) % which individual
                for ind2 = ind1:numel(inds)
                    dn1 = fullData{condition}.(inds{ind1}).(bps{bp}).(NNtests{NNtest}).dn;
                    dn2 = fullData{condition}.(inds{ind2}).(bps{bp}).(NNtests{NNtest}).dn;
                    indsName = sprintf('%svs%s',inds{ind1},inds{ind2});
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

function ksInter = ksTestInter(fullData,conditions,inds,bps,NNtests,ksInter)
% Calculates and return the individual and average ks tests inside a single
% population

% ks inter pop
% example => ksInter.condition.ind1vs1ind2.bps.NNtests = kstest2;
for NNtest = 1:numel(NNtests) % which test
    for bp = 1:numel(bps) % which brain part
        for condition1 = 1:numel(conditions) % which condition1
            for condition2 = condition1+1:numel(conditions) % which condition2
                for ind1 = 1:numel(fieldnames(fullData{condition1})) % which individual1
                    for ind2 = 1:numel(fieldnames(fullData{condition2})) % which individual2
                        % Simplify names of the 2 pops of interest
                        dn1 = fullData{condition1}.(inds{ind1}).(bps{bp}).(NNtests{NNtest}).dn;
                        dn2 = fullData{condition2}.(inds{ind2}).(bps{bp}).(NNtests{NNtest}).dn;
                        
                        % ks test per se
                        [h,p,k] = kstest2(dn1,dn2);
                        
                        % Allocate in structure
                        condsName = sprintf('%svs%s',conditions{condition1},conditions{condition2});
                        indsName = sprintf('%svs%s',inds{ind1},inds{ind2});
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

function displaySubPlots(fullData,conditions,inds,bps,NNtests,lineColors,PARAMS)
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
            for ind = 1:numel(fieldnames(fullData{condition})) % which individual
                % lot the experimental cdf
                dnExp = fullData{condition}.(inds{ind}).(bps{bp}).(NNtests{NNtest}).dn;
                %                 dn(1,:) = dn(1,:)+numel(hExp); % for tests only
                [fExp,xExp] = ecdf(dnExp); % calculate associated experimental cdf
                hExp{numel(hExp)+1} = plot(xExp,fExp,'.-','Color',lineColors(condition,:));
                % splot of the theoretical cdf => Should use dnSimu in a later version
                dnSimu = fullData{condition}.(inds{ind}).(bps{bp}).(NNtests{NNtest}).GrandCdf.mean;
                %                 dnSimu(1,:) = dnSimu(1,:)+numel(hSimu); % for tests only
                %                 [fSimu,xSimu] = ecdf(dnSimu); % calculate associated
                %                 experimental cdf => for when dnSimu = real
                fSimu = dnSimu;
                xSimu = r;
                hSimu{numel(hSimu)+1} = plot(xSimu,fSimu,'.--','Color',lineColors(condition,:));
            end
            clegs = [clegs, hExp{numel(hExp)}, hSimu{numel(hSimu)}];
            legs = [legs, conditions(condition), strcat(conditions(condition),' simu')];
        end
        legend(clegs,legs,'Location','southeast');
        legend boxoff;
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
