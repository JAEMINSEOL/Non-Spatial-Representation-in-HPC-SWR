function HistAndBar(x1,x2,CList,leg,tit)

figure('Position',[2400,300,1200,600]); 

subplot(1,3,[1 2]); hold on
histogram(x1,'Normalization','probability','facecolor',CList(1,:))
histogram(x2,'Normalization','probability','facecolor',CList(2,:))
legend(leg);
set(gca,'fontsize',12,'fontweight','b')
xlabel(tit); ylabel('proportion')

subplot(1,3,3); hold on
data = [nanmean(x1) nanmean(x2)];
err = [nanstd(x1)/sqrt(sum(~isnan(x1))) nanstd(x2)/sqrt(sum(~isnan(x2)))];

b = bar([1 2],data);
b.FaceColor = 'flat';
b.CData = CList;

    % get x positions per group
    xpos = b(1).XData + b(1).XOffset;
    % draw errorbar
    errorbar(xpos, data, err, 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
set(gca,'fontsize',12,'fontweight','b')
xticks([1 2])
xticklabels(leg)
ylabel(tit)
sgtitle(tit)