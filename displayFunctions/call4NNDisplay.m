

function call4NNDisplay()
%{

Single call function for displays of data post processing
Can call for:
- RMSE map
- exp vs CSR cdfs


load dataCombined_RMSMap
CAUTION, will load the whole file into MATLAB, easily eating 60Go of RAM

%}


%%%%%%%%%%%%%%%%%%%%%%%%% PARAMS %%%%%%%%%%%%%%%%%%%%%%%%% 
doRMSmap = 0;
doCSR = 1;
% popTestsFields = {'t2vst2', 't3vst3', 't3vst2'};
popTestsFields = {'t2vst2'};

%%%%%%%%%%%%%%%%%%%%%%%%% PARAMS %%%%%%%%%%%%%%%%%%%%%%%%% 







filePaths = uipickfiles('Prompt',...
    'Select the correct file to analyse (ex: 180201_..._dmso24h_..._yzCASE1_all..._max3cellDia_RMSMap.mat)');

for file = 1:length(filePaths)

    close all
    
    [fpath, fname, ~] = fileparts(filePaths{file});    
    
    fprintf('Processing file s%', fname);
    
    load(filePaths{file});
    
    saveFilePath = fullfile(fpath, regexprep(fname, '_RMSMap', ''));
    if doRMSmap
        RMSmap_display(dataCombined_RMSMap, saveFilePath);
    end
    if doCSR
        expVScsr_display(dataCombined_RMSMap, saveFilePath, popTestsFields);
    end
    
end

end



