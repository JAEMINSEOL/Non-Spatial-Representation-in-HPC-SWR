   units = UnitsTable_B

    figure;

    subplot(1,3,1)
    x = abs(units.RDI_LScene); y = units.RateP_LScene;
    x(isnan(y))=[]; y(isnan(y))=[];
    scatter(x,y,40,'k','filled')
    r = corrcoef(x,y);

    p = polyfit(x,y,1);
yfit = polyval(p,x);
hold on
plot(x,yfit,'r-')
text(.6,.3,['corrcoef = ' jjnum2str(r(2,1),2)],'color','r')
ylim([0 .55])
title('Left scene pair'); xlabel("Cohen's d"); ylabel('permutation p')


subplot(1,3,2)
    x = abs(units.RDI_RScene); y = units.RateP_RScene;
    x(isnan(y))=[]; y(isnan(y))=[];
    y(isnan(x))=[];x(isnan(x))=[]; 
    scatter(x,y,40,'k','filled')
    r = corrcoef(x,y);

    p = polyfit(x,y,1);
yfit = polyval(p,x);
hold on
plot(x,yfit,'r-')
text(.6,.3,['corrcoef = ' jjnum2str(r(2,1),2)],'color','r')
ylim([0 .55])
title('Right scene pair'); xlabel("Cohen's d"); ylabel('permutation p')

subplot(1,3,3)
    x = abs(units.RDI_LR); y = units.RateP_LR;
    x(isnan(y))=[]; y(isnan(y))=[];
    y(isnan(x))=[];x(isnan(x))=[]; 
    scatter(x,y,40,'k','filled')
    r = corrcoef(x,y);

    p = polyfit(x,y,1);
yfit = polyval(p,x);
hold on
plot(x,yfit,'r-')
text(.8,.3,['corrcoef = ' jjnum2str(r(2,1),2)],'color','r')
ylim([0 .55])
title('Choice')
xlabel("Cohen's d"); ylabel('permutation p')
