Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R1_300'];
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R2_300'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
RippleList_old = readtable([ROOT.Old '\RipplesTable_CA1.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
Experimenter = {'JS','LSM','SEB'};

RipplesTable = RippleList_old;

fd = dir(ROOT.Old);

RipplesTable_all = table;
thisSID_old = '';

RipplesTable_filtered = readtable([ROOT.Save '\RipplesTable_Behav_' thisRegion '_speed_filtered.xlsx']);
RipplesTable_ori = RipplesTable_filtered;
RipplesTable_f = RipplesTable_ori(RipplesTable_ori.ensemble>=3,:);
RipplesTable_JSW = readtable([ROOT.Info '\RipplesTable_JSW.xlsx']);
RipplesTable_JSW_ITI = RipplesTable_JSW(RipplesTable_JSW.area==4 & RipplesTable_JSW.ensemble>=3,:);

cl_code = {hex2rgb('712B75'),hex2rgb('56519E'),hex2rgb('1B71B7'),hex2rgb('FD8B5A')};

%%
RippleList = RipplesTable_f;
sid = 3562;
        thisRSID = [jmnum2str(RippleList.rat(sid),3) '-' jmnum2str(RippleList.session(sid),2)];
        Recording_region_TT = Recording_region({thisRSID},:);
        
        
        temp_r = decimalToBinaryVector(RippleList.TTs(sid));
        RippleTT = flip(length((temp_r))+1 - find(temp_r))';
        
        if ~strcmp(thisRSID, thisSID_old)
            TargetTT = GetTargetTT(ROOT,thisRSID,thisRegion,Params,1);
            disp([thisRSID ' EEG data loading...'])
            EEG = LoadEEGData(ROOT, thisRSID, TargetTT,Params,Params_Ripple);
            thisSID_old = thisRSID;
        end
        
        %%
        plot(EEG.TT5.Raw(RippleList.STindex(sid)-150:RippleList.EDindex(sid)+150),'color',hex2rgb('#006E4B'),'linewidth',4)
        axis off
%%
figure;
histogram((RipplesTable_f.RippleDuration)*1000,'normalization','probability','Numbins',40,'BinLimits',[0 250],'facecolor',cl_code{1})
hold on
histogram(RipplesTable_JSW_ITI.RippleDuration*1000,'normalization','probability','Numbins',40,'BinLimits',[0 250],'facecolor',cl_code{4})
title('ripple duration')
xlabel('duration(ms)')
ylabel('proportion'); ylim([0 0.15])
set(gca, 'fontsize', 15,'fontweight','b')
[h,p] = ttest2(RipplesTable_f.RippleDuration,RipplesTable_JSW_ITI.RippleDuration)
[h,p] = kstest2(RipplesTable_f.RippleDuration,RipplesTable_JSW_ITI.RippleDuration)

text(150,0.1,['Mean = ' jjnum2str(nanmean((RipplesTable_f.RippleDuration)*1000),2) 'ms'],'color',cl_code{1})
text(150,0.09,['STD = ' jjnum2str(nanstd((RipplesTable_f.RippleDuration)*1000),2) 'ms'],'color',cl_code{1})
text(150,0.075,['Mean = ' jjnum2str(nanmean(RipplesTable_JSW_ITI.RippleDuration*1000),2) 'ms'],'color',cl_code{4})
text(150,0.065,['STD = ' jjnum2str(nanstd(RipplesTable_JSW_ITI.RippleDuration*1000),2) 'ms'],'color',cl_code{4})
%%
figure;
histogram(RipplesTable_f.MeanFilteredPeak,'normalization','probability','Numbins',40,'BinLimits',[0 400],'facecolor',cl_code{1})
hold on
histogram(RipplesTable_JSW_ITI.MeanFilteredPeak,'normalization','probability','Numbins',40,'BinLimits',[0 400],'facecolor',cl_code{4})

title('LFP peak voltage')
xlabel('LFP peak voltage(150-250Hz filtered)')
ylabel('proportion')
set(gca, 'fontsize', 15,'fontweight','b')
% [h,p] = ttest2(RipplesTable_f.MeanFilteredPeak,RipplesTable_JSW_ITI.MeanFilteredPeak)
% [h,p] = kstest2(RipplesTable_f.MeanFilteredPeak,RipplesTable_JSW_ITI.MeanFilteredPeak)

text(150,0.1,['Mean = ' jjnum2str(nanmean(RipplesTable_f.MeanFilteredPeak),2) 'uV'],'color',cl_code{1})
text(150,0.09,['STD = ' jjnum2str(nanstd(RipplesTable_f.MeanFilteredPeak),2) 'uV'],'color',cl_code{1})
text(150,0.075,['Mean = ' jjnum2str(nanmean(RipplesTable_JSW_ITI.MeanFilteredPeak),2) 'uV'],'color',cl_code{4})
text(150,0.065,['STD = ' jjnum2str(nanstd(RipplesTable_JSW_ITI.MeanFilteredPeak),2) 'uV'],'color',cl_code{4})
%%
figure;
subplot(1,3,1)
hist_prefit(1e3*RipplesTable_f.RippleDuration(strcmp(RipplesTable_f.experimenter,'LSM')),cl_code{1},220,.15)
xlabel('duration (ms)')
title(['LSM ripple duration']) 

subplot(1,3,2)
hist_prefit(1e3*RipplesTable_f.RippleDuration(strcmp(RipplesTable_f.experimenter,'SEB')),cl_code{2},220,.15)
xlabel('duration (ms)')
title(['SEB ripple duration']) 

subplot(1,3,3)
hist_prefit(1e3*RipplesTable_f.RippleDuration(strcmp(RipplesTable_f.experimenter,'JS')),cl_code{3},220,.15)
xlabel('duration (ms)')
title(['JS ripple duration']) 

% [h,p] = ttest2(RipplesTable_f.RippleDuration(strcmp(RipplesTable_f.experimenter,'LSM')),...
%     RipplesTable_f.RippleDuration(strcmp(RipplesTable_f.experimenter,'SEB')))
% [p,h] = anova1(RipplesTable_f.RippleDuration*1e3,RipplesTable_f.experimenter)
% [p,h] = kruskalwallis(RipplesTable_f.RippleDuration,RipplesTable_f.experimenter)
%%
figure;
subplot(1,3,1)
hist_prefit(RipplesTable_f.MeanFilteredPeak(strcmp(RipplesTable_f.experimenter,'LSM')),cl_code{1},400,.15)
xlabel('LFP peak voltage(150-250Hz filtered)')
title(['LSM ripple peak (TT mean)']) 

subplot(1,3,2)
hist_prefit(RipplesTable_f.MeanFilteredPeak(strcmp(RipplesTable_f.experimenter,'SEB')),cl_code{2},400,.15)
xlabel('LFP peak voltage(150-250Hz filtered)')
title(['SEB ripple peak (TT mean)']) 

subplot(1,3,3)
hist_prefit(RipplesTable_f.MeanFilteredPeak(strcmp(RipplesTable_f.experimenter,'JS')),cl_code{3},400,.15)
xlabel('LFP peak voltage(150-250Hz filtered)')
title(['JS ripple peak (TT mean)']) 


% [h,p] = ttest2(RipplesTable_f.MeanFilteredPeak(strcmp(RipplesTable_f.experimenter,'SEB')),...
%     RipplesTable_f.MeanFilteredPeak(strcmp(RipplesTable_f.experimenter,'JS')))
% [p,h] = kruskalwallis(RipplesTable_f.MeanFilteredPeak,RipplesTable_f.rat)
%%
figure;
subplot(1,3,1)
hist_prefit(RipplesTable_f.StdFilteredPeak(strcmp(RipplesTable_f.experimenter,'LSM')),'r',200,.15)
xlabel('LFP peak voltage(150-250Hz filtered)')
title(['LSM ripple peak (TT std)']) 

subplot(1,3,2)
hist_prefit(RipplesTable_f.StdFilteredPeak(strcmp(RipplesTable_f.experimenter,'SEB')),'b',200)
xlabel('LFP peak voltage(150-250Hz filtered)')
ylabel('proportion')
title(['SEB ripple peak (TT std)']) 
ylim([0 0.15])
subplot(1,3,3)
hist_prefit(RipplesTable_f.StdFilteredPeak(strcmp(RipplesTable_f.experimenter,'JS')),'k',200)
xlabel('LFP peak voltage(150-250Hz filtered)')
ylabel('proportion')
title(['JS ripple peak (TT std)']) 
ylim([0 0.15])
%%
figure;
subplot(1,2,1)
boxplot((RipplesTable_f.RippleDuration)*1000)
ylim([0 250])
set(gca, 'fontsize', 15,'fontweight','b')
subplot(1,2,2)
boxplot(RipplesTable_JSW_ITI.RippleDuration*1000)
ylim([0 250])

set(gca, 'fontsize', 15,'fontweight','b')

[p,h] = ranksum(RipplesTable_f.RippleDuration,RipplesTable_JSW_ITI.RippleDuration)

quantile(RipplesTable_f.RippleDuration*1000,3)
quantile(RipplesTable_JSW_ITI.RippleDuration*1000,3)


%%
figure;
subplot(1,2,1)
boxplot(RipplesTable_f.MeanFilteredPeak)
ylim([0 350])
set(gca, 'fontsize', 15,'fontweight','b')
subplot(1,2,2)
boxplot(RipplesTable_JSW_ITI.MeanFilteredPeak)
ylim([0 350])

set(gca, 'fontsize', 15,'fontweight','b')
[h,p] = ranksum(RipplesTable_f.MeanFilteredPeak,RipplesTable_JSW_ITI.MeanFilteredPeak)

quantile(RipplesTable_f.MeanFilteredPeak,3)
quantile(RipplesTable_JSW_ITI.MeanFilteredPeak,3)


%%
figure;
subplot(1,3,1)
tempT = RipplesTable_f(strcmp(RipplesTable_f.experimenter,'LSM'),:);
boxplot(tempT.NumTT./tempT.NumAllTT)
xlabel('LSM dataset')
ylabel('% NumTT')
set(gca, 'fontsize', 15,'fontweight','b')
subplot(1,3,2)
tempT = RipplesTable_f(strcmp(RipplesTable_f.experimenter,'SEB'),:);
boxplot(tempT.NumTT./tempT.NumAllTT)
xlabel('Delcasso dataset')
set(gca, 'fontsize', 15,'fontweight','b')
subplot(1,3,3)
tempT = RipplesTable_f(strcmp(RipplesTable_f.experimenter,'JS'),:);
boxplot(tempT.NumTT./tempT.NumAllTT)
xlabel('Jhoseph dataset')
set(gca, 'fontsize', 15,'fontweight','b')



%%
% Add Ripple quality

RipplesTable_f.SPKperSEC = RipplesTable_f.spike ./ RipplesTable_f.RippleDuration;
RipplesTable_f.SPKperENS = RipplesTable_f.spike ./ RipplesTable_f.ensemble;
RipplesTable_JSW_ITI.SPKperSEC = RipplesTable_JSW_ITI.spike ./ RipplesTable_JSW_ITI.RippleDuration;
RipplesTable_JSW_ITI.SPKperENS = RipplesTable_JSW_ITI.spike ./ RipplesTable_JSW_ITI.ensemble;
%%
figure;
histogram(RipplesTable_f.SPKperSEC,'normalization','probability','Numbins',40,'BinLimits',[0 800],'facecolor',cl_code{1})
hold on
histogram(RipplesTable_JSW_ITI.SPKperSEC,'normalization','probability','Numbins',40,'BinLimits',[0 800],'facecolor',cl_code{4})
title('spike/sec distribution')
xlabel('spike / sec')
ylabel('proportion')
set(gca, 'fontsize', 15,'fontweight','b')
[h,p]=ttest2(RipplesTable_f.SPKperSEC,RipplesTable_JSW_ITI.SPKperSEC)
%%
figure;
subplot(1,2,1)
boxplot(RipplesTable_f.SPKperSEC)
ylim([0 800])
set(gca, 'fontsize', 15,'fontweight','b')
subplot(1,2,2)
boxplot(RipplesTable_JSW_ITI.SPKperSEC)
ylim([0 800])

set(gca, 'fontsize', 15,'fontweight','b')
[h,p] = ranksum(RipplesTable_f.SPKperSEC,RipplesTable_JSW_ITI.SPKperSEC)

quantile(RipplesTable_f.SPKperSEC,3)
quantile(RipplesTable_JSW_ITI.SPKperSEC,3)

%%
figure;
subplot(2,2,1)
hist_prefit(RipplesTable_f.SPKperSEC(strcmp(RipplesTable_f.experimenter,'LSM')),cl_code{1},800,.22)
xlabel('spike / sec')
title(['LSM spike / sec']) 

subplot(2,2,2)
hist_prefit(RipplesTable_f.SPKperSEC(strcmp(RipplesTable_f.experimenter,'SEB')),cl_code{2},800,.22)
xlabel('spike / sec')
title(['Delcasso spike / sec']) 

subplot(2,2,3)
hist_prefit(RipplesTable_f.SPKperSEC(strcmp(RipplesTable_f.experimenter,'JS')),cl_code{3},800,.22)
xlabel('spike / sec')
title(['Jhoseph spike / sec']) 

subplot(2,2,4)
hist_prefit(RipplesTable_JSW_ITI.SPKperSEC,cl_code{4},800,.22)
xlabel('spike / sec')
title(['JSW spike / sec']) 

%%
figure;
histogram(RipplesTable_f.SPKperENS,'normalization','probability','Numbins',40,'BinLimits',[1 6],'facecolor',cl_code{1})
hold on
histogram(RipplesTable_JSW_ITI.SPKperENS,'normalization','probability','Numbins',40,'BinLimits',[1 6],'facecolor',cl_code{4})
title('spike / cell distribution')
xlabel('spike / cell')
ylabel('proportion')
set(gca, 'fontsize', 15,'fontweight','b')
[h,p]=ttest2(RipplesTable_f.SPKperENS,RipplesTable_JSW_ITI.SPKperENS)
%%
figure;
subplot(1,2,1)
boxplot(RipplesTable_f.SPKperENS)
ylim([0 6])
set(gca, 'fontsize', 15,'fontweight','b')
subplot(1,2,2)
boxplot(RipplesTable_JSW_ITI.SPKperENS)
ylim([0 6])

set(gca, 'fontsize', 15,'fontweight','b')
[h,p] = ranksum(RipplesTable_f.SPKperENS,RipplesTable_JSW_ITI.SPKperENS)

quantile(RipplesTable_f.SPKperENS,3)
quantile(RipplesTable_JSW_ITI.SPKperENS,3)

%%
figure;
subplot(2,2,1)
hist_prefit(RipplesTable_f.SPKperENS(strcmp(RipplesTable_f.experimenter,'LSM')),cl_code{1},6,.31)
xlabel('spike / cell')
title(['LSM spike / cell']) 

subplot(2,2,2)
hist_prefit(RipplesTable_f.SPKperENS(strcmp(RipplesTable_f.experimenter,'SEB')),cl_code{2},6,.31)
xlabel('spike / cell')
title(['Delcasso spike / cell']) 

subplot(2,2,3)
hist_prefit(RipplesTable_f.SPKperENS(strcmp(RipplesTable_f.experimenter,'JS')),cl_code{3},6,.31)
xlabel('spike / cell')
title(['Jhoseph spike / cell']) 

subplot(2,2,4)
hist_prefit(RipplesTable_JSW_ITI.SPKperENS,cl_code{4},6,.31)
xlabel('spike / cell')
title(['JSW spike / cell']) 



%%
function histfit_n(data,cl)
hh = histfit(data,40);
hh(1).FaceAlpha = 0; hh(1).EdgeAlpha = 0; hh(2).Color = 'w';
ybar = hh(1).YData;          
plot(hh(2).XData, hh(2).YData/sum(ybar), cl, 'LineWidth',2)

end

function hist_prefit(data,cl,x_lim,y_lim)
hh = histogram(data,'normalization','probability','Numbins',40,'BinLimits',[1 x_lim]);
hh.FaceColor = cl;
text(x_lim*0.7,0.06,['N = ' jjnum2str(length(data),2)],'color',cl)
text(x_lim*0.7,0.05,['Mean = ' jjnum2str(nanmean(data),2)],'color',cl)
text(x_lim*0.7,0.04,['Std = ' jjnum2str(nanstd(data),2)],'color',cl)
set(gca, 'fontsize', 15,'fontweight','b')
ylabel('proportion')
ylim([0 y_lim])
end
