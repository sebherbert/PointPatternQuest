
function errorRMSmap()
% displays the RMSE map 

% Calculate hessian matrix and second order partial derivatives at the
% point of interest
% only works for t3vst2 for the moment, could be expanded

% Pick all the folders containing files of interest
listFolders = uipickfiles('Prompt','Select folders to process');

% Enter the proper regexp for file import
prompt = {'Enter RMS map regexp:','Enter fitted model regexp:'};
dlg_title = 'Input';
num_lines = 1;
defaultans = {'*_dm_t3vst2_RMSMap*','*_dm_fit*'};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

fileRMSMapregexp = answer{1};
fileFitregexp = answer{2};

% Loop through the different folders
for folder = 1:numel(listFolders)
    
    fprintf('Processing folder: %s\n', listFolders{folder});
    
    close all
    % Go to the folder of interest
    cd(listFolders{folder})
    
    % Check and load the fit file
    fitFile = {};
    tempFitFile = dir (fileFitregexp);
    for file = 1:numel(tempFitFile)
        if ~contains(tempFitFile(file).name,'_delta')
            fitFile{length(fitFile)+1} = tempFitFile(file);
        end
    end
    if numel(fitFile)~=1
        fprintf('Error too many fits in folder: %s\n',listFolders{folder});
        continue
    else
        fitData = load(fitFile{1}.name);
        fitField = fieldnames(fitData);
    end
    
    % Check and load the RMS file
    RMSfile = dir (fileRMSMapregexp);
    if numel(RMSfile)~=1
        fprintf('Error too many RMSE files in folder: %s\n',listFolders{folder});
        continue
    else
        RMSdata = load(RMSfile.name);
        RMSField = fieldnames(RMSdata);
    end
    
    assoDer = evaluateError2ndDer(RMSdata.(RMSField{1}),fitData.(fitField{1}),...
        RMSfile,fitFile);
    fitData.(fitField{1}).t3vst2.assoDer = assoDer;
    if fitField{1} == 'dataCombined'
        dataCombined = fitData.(fitField{1});
    else
        fprintf('Error in folder: %s. Expecting a dataCombined field\n',listFolders{folder});
        continue
    end
    [~, fileName, ~] = fileparts(fitFile{1}.name);
    save(strcat(fileName,'_deltaSmoothImguidedfilter.mat'),'dataCombined');

    clear RMSdata fitData
end

end

function assoDer = evaluateError2ndDer(RMSdata,fitData,RMSfile,fitfile)
%% format finalRMS as a 2D map
% extract and concatenate fields of interest

% get list of models
fieldList = fieldnames(RMSdata);

% get table of results
allModels = table;
for field = 1:numel(fieldList) % for each field
    % Check if field is a model
    if contains(fieldList{field},'Model')
        allModels = vertcat(allModels, ...
            table(RMSdata.(fieldList{field}).effectRange, ...
            RMSdata.(fieldList{field}).effectStrength, ...
            RMSdata.(fieldList{field}).t3vst2.medRMS));
    end
end
variableNames = {'Range','Strength','finalRMS'};
allModels.Properties.VariableNames = variableNames;

allR = unique(allModels.Range);
allS = unique(allModels.Strength);

% Reshape
mapRMS = reshape(allModels.finalRMS,[numel(allS),numel(allR)]);
allRmeshForm = reshape(allModels.Range,[numel(allS),numel(allR)]);
allSmeshForm = reshape(allModels.Strength,[numel(allS),numel(allR)]);


%% improve smoothness
% mapRMS = imgaussfilt(mapRMS,0.5);
% mapRMS = imguidedfilter(mapRMS);

%% Interpolate along the Rf and Sf at exact fitted values

Rf = fitData.t3vst2.fit.Range;
Sf = fitData.t3vst2.fit.Strength;

NInter = 5; % Number of steps used for interpolation

% Calculation of R and S steps (S steps are not linear!)
% Epsilon is estimated for 1/Nth of a step around the fitted value
epsilonStep = 2;
Rstep = (allR(2)-allR(1))/epsilonStep;

temp = allS - Sf;
Sstep = (allS(numel(temp(temp<0))+1) - allS(numel(temp(temp<0))))/epsilonStep;

% interpolation along Rf
Rinter = Rf-NInter*Rstep : Rstep : Rf+NInter*Rstep; % NInter steps around the fitted value
RMSf_interR = interp2(allR, allS, mapRMS, Rinter, Sf, 'spline');

% interpolation along Sf
Sinter = Sf-NInter*Sstep : Sstep : Sf+NInter*Sstep; % NInter steps around the fitted value
RMSf_interS = interp2(allR, allS, mapRMS, Rf, Sinter, 'spline');

% interpolation at Rf and Sf
RMSf = interp2(allR, allS, mapRMS, Rf, Sf);

%% Calculate seconde order gradient along R and S.
% By interpolating specific additionnal positions
gradR = gradient(RMSf_interR, Rstep);
gradRR = gradient(gradR, Rstep);
errorR = real(sqrt(1./gradRR));
errorR(errorR==0) = nan;
errorRf = errorR(ceil(numel(errorR)/2));

gradS = gradient(RMSf_interS, Sstep);
gradSS = gradient(gradS, Sstep);
errorS = real(sqrt(1./gradSS));
errorS(errorS==0) = nan;
errorSf = errorS(ceil(numel(errorS)/2));

% without interpolation at individual line level, although interpolation is still possible! Just
% don't specify the positions of the new positions and remember that your
% points of measure are already FUCKING SEMI-LOG!!!!

interFact = 2;
newMapRMS = interp2(mapRMS,interFact,'spline'); % multiply each dimension by 4
newAllR = linspace(allR(1),allR(end),(numel(allR)-1)*2^interFact+1);
newAllS = exp(linspace(log(allS(1)),log(allS(end)),(numel(allS)-1)*2^interFact+1));

% full map of the gradients
newRstep = newAllR(2)-newAllR(1);
newSstep = log(newAllS(2))-log(newAllS(1));
gr = gradient(double(newMapRMS), newRstep);
gs = gradient(double(newMapRMS), newSstep);
grr = gradient(gr, newRstep);
grs = gradient(gr, newSstep);
gsr = gradient(gs, newRstep);
gss = gradient(gs, newSstep);

% => does not take into account step size (old version)
[gr, gs] = gradient(double(newMapRMS), newSstep, newRstep);
[grr, grs] = gradient(gr, newSstep, newRstep);
[gsr, gss] = gradient(gs, newSstep, newRstep);

RMS2nd_Rf = interp2(newAllR, log(newAllS), grr, Rf, log(Sf), 'spline');
RMS2nd_Sf = interp2(newAllR, log(newAllS), gss, Rf, log(Sf), 'spline');

RMS_Rferror = sqrt(1/RMS2nd_Rf);
RMS_Sferror = sqrt(1/RMS2nd_Sf);


%% Save datafile
assoDer = {};
assoDer.RMSfile = RMSfile;
assoDer.fitfile = fitfile;
assoDer.allModels = allModels;
assoDer.Rf = Rf;
assoDer.Sf = Sf;
% assoDer.NInter = NInter;
% assoDer.epsilonStep = epsilonStep;
% assoDer.errorR = errorR;
% assoDer.errorS = errorS;
% assoDer.errorRf = errorRf;
% assoDer.errorSf = errorSf;
assoDer.interFact = interFact;
assoDer.maps.gr = gr;
assoDer.maps.gs = gs;
assoDer.maps.grr = grr;
assoDer.maps.grs = grs;
assoDer.maps.gss = gss;
assoDer.maps.gsr = gsr;
assoDer.fast.RMS2nd_Rf = RMS2nd_Rf;
assoDer.fast.RMS2nd_Sf = RMS2nd_Sf;
assoDer.fast.RMS_Rferror = RMS_Rferror;
assoDer.fast.RMS_Sferror = RMS_Sferror;
save('associatedError','assoDer');

%% Displays
% display surface and fitted parameters
figure; 
surf(allR, allS, mapRMS, 'EdgeColor', 'None', 'FaceColor', 'interp');
% plot3(allModels.Range, log(allModels.Strength), allModels.finalRMS)
% display the fitted position
hold on
scatter3(Rf,Sf,RMSf,'filled','MarkerFaceColor',[217 83 25]/255)
xlabel('Range'); ylabel('Strength'); zlabel('RMSE');
set(gca, 'YScale', 'log');
saveas(gcf,'RMSmap_plus_fit.png');
saveas(gcf,'RMSmap_plus_fit.fig');

% % display full map
% figure
% hold on
% surf(newAllR, newAllS, newMapRMS, 'EdgeColor', 'None', 'FaceColor', 'interp');
% surf(newAllR, newAllS, gr, 'EdgeColor', 'None', 'FaceColor', 'interp');
% surf(newAllR, newAllS, grr, 'EdgeColor', 'None', 'FaceColor', 'interp');
% scatter3(Rf,Sf,RMSf,'filled','MarkerFaceColor',[217 83 25]/255)
% scatter3(Rf,Sf,RMS2nd_Rf,'filled','MarkerFaceColor',[217 83 25]/255)
% 
% xlabel('Range'); ylabel('Strength');
% set(gca, 'YScale', 'log');
% legend({'RMS map' 'Gradient along R' 'Gradient 2nd along R'})
% saveas(gcf,'RMSmap_plus_gradsR.png');
% saveas(gcf,'RMSmap_plus_gradsR.fig');

% % display full map
% figure
% hold on
% surf(newAllR, newAllS, newMapRMS, 'EdgeColor', 'None', 'FaceColor', 'interp');
% surf(newAllR, newAllS, gs, 'EdgeColor', 'None', 'FaceColor', 'interp');
% surf(newAllR, newAllS, gss, 'EdgeColor', 'None', 'FaceColor', 'interp');
% scatter3(Rf,Sf,RMSf,'filled','MarkerFaceColor',[217 83 25]/255)
% scatter3(Rf,Sf,RMS2nd_Sf,'filled','MarkerFaceColor',[217 83 25]/255)
% xlabel('Range'); ylabel('Strength');
% set(gca, 'YScale', 'log');
% legend({'RMS map' 'Gradient along S' 'Gradient 2nd along S'})
% saveas(gcf,'RMSmap_plus_gradsS.png');
% saveas(gcf,'RMSmap_plus_gradsS.fig');

% % Example of a fit in S
% figure
% hold on
% plot(RMSf_interS)
% plot(gradS)
% plot(gradSS)
% yyaxis right
% plot(errorS);
% legend({'RMS interpolation S' 'Gradient along S' 'Second order gradient' 'Local measurement error {\surd} (1/f'''')'})
% saveas(gcf,'1D_example.png');
% saveas(gcf,'1D_example.fig');






end



