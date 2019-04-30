
function expVScsr_display(dataCombined_RMSMap, saveFilePath)
%% Display a CDF of the experimental data vs Complete Spatial Randomness
% version v0p1

% load ...dataCombined_RMSMap => Contains the CSR among other models
% All the 
% CAUTION, will load the whole file into MATLAB, easily eating 60Go of RAM

% Reassociate PARAMS
PARAMS = dataCombined_RMSMap.PARAMS;


% get list of models
fieldList = fieldnames(dataCombined_RMSMap);
popTestsFields = {'t2vst2', 't3vst3', 't3vst2'};

% get table of results
allCSRs = {};
for field = 1:numel(fieldList) % for each field
    % Check if field is a NONE model
    if contains(fieldList{field},'Model_TN') % Should contain all the NONE (effect) models

        for popTestsNum = 1:length(popTestsFields) % For each pop cross test
            
            popCross = popTestsFields{popTestsNum};
            
            if ~isfield(allCSRs, popCross)
                % Initialize dnSimu, dnExp, simuCDFs and expCDFs fields
                allCSRs.(popCross).dnSimu = [];
                allCSRs.(popCross).dnExp = [];
                allCSRs.(popCross).simuCDFs = [];
                allCSRs.(popCross).expCDFs = [];
                
                % Save dnExp (Experimental distribution) at the first encounter only
                allCSRs.(popCross).dnExp = ...
                    dataCombined_RMSMap.(fieldList{field}).(popCross).dnExp;
                allCSRs.(popCross).expCDFs = ...
                    dataCombined_RMSMap.(fieldList{field}).(popCross).expCDFs;
                
                % Save some results specific to this popCross
                allCSRs.(popCross).PARAMS = ...
                    dataCombined_RMSMap.(fieldList{field}).(popCross).PARAMS;
            
            end
            
            allCSRs.(popCross).dnSimu = ...
                horzcat(dataCombined_RMSMap.(fieldList{field}).(popCross).dnSimu, allCSRs.(popCross).dnSimu);
        end
    end
end

% recalculate CDFs
for popTestsNum = 1:length(popTestsFields) % For each pop cross test
    popCross = popTestsFields{popTestsNum};
    allCSRs.(popCross).simuCDFs = formatCdfsSimu(allCSRs.(popCross).dnSimu, PARAMS);
end

% display CDFs
for popTestsNum = 1:length(popTestsFields) % For each pop cross test
    popCross = popTestsFields{popTestsNum};
    figure
    displayNNCDFs(allCSRs.(popCross).expCDFs, allCSRs.(popCross).simuCDFs, ...
        allCSRs.(popCross).PARAMS, allCSRs.(popCross).PARAMS.useRMSMaxDist);
    saveas(gcf,sprintf('%s_%s_exp_VS_CSR', saveFilePath, popCross));
    saveas(gcf,sprintf('%s_%s_exp_VS_CSR.png', saveFilePath, popCross));
end

end


function simuCDFs = formatCdfsSimu(dnSimu, PARAMS)
% New version using a "simpler" and more manual percentile approach (still kept 
% Greenwood for the experimental cdf. (same process as original WS's)

% Calculate the CDFs of the simulated populations
for simu = 1:size(dnSimu,2) % each simulation is treated individually
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




