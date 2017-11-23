

function GrandAll = displaySimuPPQ(xScale, Grand, NN, NumPermut, figTitle, figSavePath)
% display the figure bearing the Nearest Neighbours histogram and the
% random distribution simulations

GrandAll = {};
for i=1:size(xScale,2)
    GrandAll.mean(i)=mean(Grand(:,i));
    GrandAll.std(i)=std(Grand(:,i));
    GrandAll.iqr5(i)=prctile(Grand(:,i),5);
    GrandAll.iqr95(i)=prctile(Grand(:,i),95);
    GrandAll.iqr1(i)=prctile(Grand(:,i),1);
    GrandAll.iqr99(i)=prctile(Grand(:,i),99);
end


figure
ylabel('Cumulative cell frequency');
xlabel('distance to nearest neighbor (in cell diameter)');

hold on
plot(xScale,GrandAll.mean,'-k','linewidth',4);
plot(xScale,GrandAll.iqr5,'-','linewidth',3,'color',[0.6 0.6 0.6]);
plot(xScale,GrandAll.iqr95,'-','linewidth',3,'color',[0.6 0.6 0.6]);
plot(xScale,GrandAll.iqr1,'-','linewidth',1,'color',[0.6 0.6 0.6]);
plot(xScale,GrandAll.iqr99,'-','linewidth',1,'color',[0.6 0.6 0.6]);

plot(xScale,NN,'-r','linewidth',4);
text(0.5,0.95,'95% and 99% intervals')
text(0.5,0.9,[num2str(NumPermut),' random perm.'])

% force the axes
if 2*ceil(xScale(min(find(NN>0.9))))<14
    axis([0 2*ceil(xScale(min(find(NN>0.9)))) 0 1]);
else
    axis([0 14 0 1]);
end

title(figTitle,'FontSize',10)

saveas(gcf,figSavePath);
% saveas(gcf,figSavePath,'tiffn');
saveas(gcf,figSavePath,'png');

end