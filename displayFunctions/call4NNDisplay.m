

function call4NNDisplay()
%{

Single call function for displays of data post processing
Can call for:
- RMSE map
- exp vs CSR cdfs


load dataCombined_RMSMap
CAUTION, will load the whole file into MATLAB, easily eating 60Go of RAM

%}

filePaths = uipickfiles('Prompt',...
    'Select the correct file to analyse (ex: 180201_..._dmso24h_..._yzCASE1_all..._max3cellDia_RMSMap.mat)');

close all

for file = 1:length(filePaths)
    
    [fpath, fname, ~] = fileparts(filePaths{file});
    load(filePaths{file});
    
    saveFilePath = fullfile(fpath, regexprep(fname, '_RMSMap', ''));
    RMSmap_display(dataCombined_RMSMap, saveFilePath);
    expVScsr_display(dataCombined_RMSMap, saveFilePath);
    
end

end



