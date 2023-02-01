function fig = popul_FRMap(FRMap,yLabel,xTicks,xLabel,hand)


fig=imagesc(FRMap);

axis off
for i=1:length(xTicks)
    text(xTicks(i),size(FRMap,1)+1,xLabel{i},'color','k')
end
clim([0 1])
if hand
for i=1:size(FRMap,1)
    text(-8,i,yLabel{size(FRMap,1)-i+1},'color','k')
end
end