

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
% => Objective is to create 3 structures with each 3 samples, then 3 brain
% parts and then 3 tests
% pop1Files = ;

load('sox2_C_subdiv_L_corrected_nodb_noDl_xyzCASE1_allSample.mat'); % temporary

dataYoung.ind1.all = dataCombined;
dataYoung.ind2.all = dataCombined;
dataYoung.ind3.all = dataCombined;
dataOld.ind1.all = dataCombined;
dataOld.ind2.all = dataCombined;
dataOld.ind3.all = dataCombined;
dataDrug.ind1.all = dataCombined;
dataDrug.ind2.all = dataCombined;
dataDrug.ind3.all = dataCombined;

dataYoung.ind1.da = dataCombined;
dataYoung.ind2.da = dataCombined;
dataYoung.ind3.da = dataCombined;
dataOld.ind1.da = dataCombined;
dataOld.ind2.da = dataCombined;
dataOld.ind3.da = dataCombined;
dataDrug.ind1.da = dataCombined;
dataDrug.ind2.da = dataCombined;
dataDrug.ind3.da = dataCombined;

dataYoung.ind1.dl = dataCombined;
dataYoung.ind2.dl = dataCombined;
dataYoung.ind3.dl = dataCombined;
dataOld.ind1.dl = dataCombined;
dataOld.ind2.dl = dataCombined;
dataOld.ind3.dl = dataCombined;
dataDrug.ind1.dl = dataCombined;
dataDrug.ind2.dl = dataCombined;
dataDrug.ind3.dl = dataCombined;

%% Initialize the parameters
inds = {'ind1','ind2','ind3'}; % names of the individual fishes
bps = {'all','da','dl'}; % parts of the brain
NNtests = {'t2vst2','t3vst3','t3vst2'}; % tested cell types

% %% Display inter and intra (1 color each) subplot 3(parts)x3tests
% conditions = {'Young','Old'}; % Conditions to be tested
% lineColors = lines(numel(conditions));
% 
% displaySubPlots({dataYoung,dataOld},conditions,inds,bps,NNtests,lineColors);


%% ks tests
ks intra population
conditions = {'Young','Old','Drug'}; % Conditions to be tested
ksIntra = ksTestIntra({dataYoung,dataOld,dataDrug},conditions,inds,bps,NNtests);
% 
% % ks inter population
% conditions = {'Young','Old'}; % Conditions to be tested
% ksInter = ksTestInter({dataYoung,dataOld},conditions,inds,bps,NNtests);

% Reshape for reading in 3x3 tables cd excel sheets (name rows and lines)



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
                dn1 = fullData{condition}.(inds{ind1}).(bps{bp}).(NNtests{NNtest}).dn;
                for ind2 = ind1:numel(inds)
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

% function ksTestInter(fullData,conditions,inds,bps,NNtests)
% % Calculates and return the individual and average ks tests inside a single
% % population
% 
% % ks inter pop
% ksInter = {};
% % example => ksInter.condition.ind1vs1ind2.bps.NNtests = kstest2;
% for NNtest = 1:numel(NNtests) % which test
%     for bp = 1:numel(bps) % which brain part
%         for condition = 1:numel(conditions) % which condition
%             for ind1 = 1:numel(fieldnames(fullData{condition})) % which individual
% 
% 
% end
% 
% end

function displaySubPlots(fullData,conditions,inds,bps,NNtests,lineColors)
% Display inter and intra (1 color each) subplot 3(parts of brain)x3 tests
close all

figure

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
                dn = fullData{condition}.(inds{ind}).(bps{bp}).(NNtests{NNtest}).dn;
                %                 dn(1,:) = dn(1,:)+numel(hExp); % for tests only
                [fExp,xExp] = ecdf(dn); % calculate associated experimental cdf
                hExp{numel(hExp)+1} = plot(xExp,fExp,'.-','Color',lineColors(condition,:));
                % splot of the theoretical cdf => Should use dnSimu in a later version
                dnSimu = fullData{condition}.(inds{ind}).(bps{bp}).(NNtests{NNtest}).GrandAll.mean;
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
            
        
% [f,x] = ecdf(dn);
% hold on
% plot(x,f,'.-');
% plot(r,G,'.-');
% % plot(r,GrandAll.mean);

end
