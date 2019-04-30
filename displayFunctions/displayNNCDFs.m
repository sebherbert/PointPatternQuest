
function displayNNCDFs(expCDFs, simuCDFs, PARAMS, dispMaxCdf)

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
    h(7) = line([simuCDFs.x(xPos), simuCDFs.x(xPos)], ...
        [0 1], ylim, 'Color', [0 127 0]/255, 'LineStyle', ':', 'LineWidth', 2);
    % Adapt the legend
    legs = [legs sprintf('maxCDF = %0.1fum(<=>%dcellDia)',...
        PARAMS.maxDistFactor*PARAMS.cellDiameter, PARAMS.maxDistFactor)];
end

% force the axes
axis(PARAMS.axis);

% display legend
legend(h,legs,'Location','southeast');
legend boxoff;

end