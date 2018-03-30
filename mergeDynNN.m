

function outNNdispersion = mergeDynNN(PARAMS, inNNdistances)
%{
 This function combines and analyse the nearest neighbour distances
obtained by comparing the dispersion of cell types at different time points
(simulation and experiment). It takes in input the NNdistances output
variable from the extractDynNN function (one RS pair only) and outputs
multiple statistical evaluation of the distribution combined over all the
time points.

INPUT
PARAMs - See pointDynPatternMain
inNNdistances - Structure. Each field contains the tpX vs tpY NN distances.
   See function extractDynNN() variable outNNdistances

OUTPUT
outNNdispersion - Structure. Each contain a certain Δt further containing all
   relevant statistical analaysis

%}

% create output function
outNNdispersion = {};

% Recover all crossed timepoints 
crossTpList = fieldnames(inNNdistances);

% For each crossed timepoint
for crossTp = 1:numel(crossTpList)
    
    % Extract Deltat => WARNING ONLY WORKS WITH LESS THAN 10 tps => else will
    % require a more advanced parsing
    if PARAMS.movie.maxTp < 9    
        dts = crossTpList{crossTp}(strfind(crossTpList{crossTp},'dt')+2);
        dtStart = str2num(dts(1));
        dtStop = str2num(dts(2));
        deltaT = sprintf('deltaT%d',dtStop-dtStart);
    else
        msg = fprintf('Error: Too many time points. Adapt code for multiple digit frames.\n');
        error(msg);
    end
    
    if isfield(outNNdispersion, deltaT)
        outNNdispersion.(deltaT).allDistances = vertcat(outNNdispersion.(deltaT).allDistances, ...
            inNNdistances.(crossTpList{crossTp}));
    else
        outNNdispersion.(deltaT).allDistances = inNNdistances.(crossTpList{crossTp});
    end
end

end



