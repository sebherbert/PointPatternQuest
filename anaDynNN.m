

function outNNana = anaDynNN(PARAMS, expDists,simuDists)

%{ 
anaDynNN: Compares the experimental and simulated populations, the 2 inputs
share the same structure of N fields for N DeltaN and a subfield of
allDistances

INPUT
expDists - See outNNdispersion variable in the mergeDynNN function.
   Containes the merged distances for all timepoints where DeltaT is the
   desired one.
simuDists - See outNNdispersion variable in the mergeDynNN function.
   Containes the merged distances for all timepoints where DeltaT is the
   desired one. For one RS simulation only

OUTPUT
outNNana - Structure containing all the new analyses sorted by type (RMSE,
   other) and then by deltaT

%}

% % for dev only
% expDists = NNdispersion.exp;
% simuDists = NNdispersion.simu.(pairString);

% Initialize structures
outNNana = {};
RMSEstruct = {};

DTfields = fieldnames(simuDists);
for DTfield = 1:numel(DTfields) % for each DeltaT
    
    %     % if using filters before analysis => try with pchip (beware of first false line output...)
    %     [foo, bar] = ecdf(reshape(inNNdispersionSimu.deltaT4.allDistances,[],1));
    %     tempFormSimu(:,1) = foo;
    %     tempFormSimu(:,2) = bar;
    %     temp = pchip(tempFormSimu(:,2),tempFormSimu(:,1),0:150);
    %     plot(0:150,temp)
    %
    
    %% RMSE calculations    
    if PARAMS.anaGlobal.doRMSE
        % Calculate simulated and experimental cdfs
        [RMSEstruct.expY, RMSEstruct.expX] = ecdf(reshape(expDists.(DTfields{DTfield}).allDistances,[],1));
        RMSEstruct.expX(1) = 0;
        [RMSEstruct.simuY, RMSEstruct.simuX] = ecdf(reshape(simuDists.(DTfields{DTfield}).allDistances,[],1));
        RMSEstruct.simuX(1) = 0;

        % Interpolate the simulated data against the experimental ones, could do
        % the opposite or use a single array of positions using the same function
        % (better than spline in that particular case)
        RMSEstruct.interpSimuY = pchip(RMSEstruct.simuX,RMSEstruct.simuY,RMSEstruct.expX);
        
        RMSEstruct.diffCdf = RMSEstruct.expY - RMSEstruct.interpSimuY;
        
        RMSEstruct.RMSE = rms(RMSEstruct.diffCdf);
        
        outNNana.RMSE.(DTfields{DTfield}) = RMSEstruct;
        
        % To display for construction
        %         figure
        %         plot(RMSEstruct.expX,RMSEstruct.expY)
        %         hold on
        %         plot(RMSEstruct.simuX,RMSEstruct.simuY)
        %         plot(RMSEstruct.expX,RMSEstruct.smoothedSimuY, 'o')
    end    
end







end

























