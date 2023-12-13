%% UnitSpec_compare
Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed ''];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R2'];
ROOT.Rip4 = [ROOT.Save '\ripples_mat\R4_SUB_refCA1'];
ROOT.Rip5 = [ROOT.Save '\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Save '\ripples_mat\ProfilingSheet\R25_ca1'];
ROOT.Unit1 = [ROOT.Save '\units_mat\U1'];
ROOT.Units = [ROOT.Save '\units_mat\U2'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];


CList = [ [207 8 23]/255;[23 84 181]/255];


RegionList = {'SUB','CA1'};

RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_SUB_forAnalysis' '.xlsx']);
ReactTable.SUB = readtable([ROOT.Save '\ReactTable_SUB_SUB.xlsx']);
UnitsTable.SUB = readtable([ROOT.Save '\UnitsTable_SUB_forAnalysis_TP.xlsx']);
UnitsTable_field.SUB = readtable([ROOT.Save '\UnitsTable_SUB_field_forAnalysis.xlsx']);
% UnitPair.SUB = readtable([ROOT.Save '\UnitPair_SUB.xlsx']);
% UnitPair_field.SUB = readtable([ROOT.Save '\UnitPair_SUB_field.xlsx']);

RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx']);
ReactTable.CA1 = readtable([ROOT.Save '\ReactTable_CA1_CA1.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_forAnalysis_TP.xlsx']);
UnitsTable_field.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx']);
% UnitPair.CA1= readtable([ROOT.Save '\UnitPair_CA1.xlsx']);
% UnitPair_field.CA1 = readtable([ROOT.Save '\UnitPair_CA1_field.xlsx']);

Un_SUB = UnitsTable.SUB;
Un_CA1= UnitsTable.CA1;
UnF_SUB = UnitsTable_field.SUB;
UnF_CA1= UnitsTable_field.CA1;

Un_SUB(Un_SUB.rat==232 & Un_SUB.session==4,:)=[];
UnF_SUB(UnF_SUB.rat==232 & UnF_SUB.session==4,:)=[];
Un_CA1(Un_CA1.rat==232 & Un_CA1.session==4,:)=[];
UnF_CA1(UnF_CA1.rat==232 & UnF_CA1.session==4,:)=[];

%%
UnF = UnF_CA1; Un = Un_CA1;
for uid = 1:size(Un,1)
    thisID = Un.ID{uid};
    idx = strncmp(thisID,UnF.ID,12); thisF = UnF(idx,:);
    [~,t] = max(abs(thisF.RDI_LScene)); Un.RDI_LScene_field(uid) = thisF.RDI_LScene(t);
    [~,t] = max(abs(thisF.RDI_RScene)); Un.RDI_RScene_field(uid) = thisF.RDI_RScene(t);
    [~,t] = max(abs(thisF.RDI_LR)); Un.RDI_LR_field(uid) = thisF.RDI_LR(t);
end
   
UnF_CA1 = UnF; Un_CA1 = Un;

UnF = UnF_SUB; Un = Un_SUB;
for uid = 1:size(Un,1)
    thisID = Un.ID{uid};
    idx = strncmp(thisID,UnF.ID,12); thisF = UnF(idx,:);
    [~,t] = max(abs(thisF.RDI_LScene)); Un.RDI_LScene_field(uid) = thisF.RDI_LScene(t);
    [~,t] = max(abs(thisF.RDI_RScene)); Un.RDI_RScene_field(uid) = thisF.RDI_RScene(t);
    [~,t] = max(abs(thisF.RDI_LR)); Un.RDI_LR_field(uid) = thisF.RDI_LR(t);
end
   
UnF_SUB = UnF; Un_SUB = Un;

SS=Un_SUB(Un_SUB.NumField==1,:); SM=Un_SUB(Un_SUB.NumField>1,:);
CS=Un_CA1(Un_CA1.NumField==1,:); CM=Un_CA1(Un_CA1.NumField>1,:);

dat1 = [sum(nanmax([abs(SS.RDI_LScene_field) abs(SS.RDI_RScene_field) abs(SS.RDI_LR_field)],[],2)>=0.1)/size(SS,1) ...
    sum(nanmax([abs(SM.RDI_LScene_field) abs(SM.RDI_RScene_field) abs(SM.RDI_LR_field)],[],2)>=0.1)/size(SM,1)];

dat2 = [sum(nanmax([abs(CS.RDI_LScene_field) abs(CS.RDI_RScene_field) abs(CS.RDI_LR_field)],[],2)>=0.1)/size(CS,1) ...
    sum(nanmax([abs(CM.RDI_LScene_field) abs(CM.RDI_RScene_field) abs(CM.RDI_LR_field)],[],2)>=0.1)/size(CM,1)];

varList = {'RDI_LScene_field','RDI_RScene_field','RDI_LR_field'};



var1='RDI_LScene_field';

dataL = [[sum(abs(SS.(var1))>=0.1)/size(SS,1) sum(abs(SM.(var1))>=0.1)/size(SM,1)],...
    [sum(abs(CS.(var1))>=0.1)/size(CS,1)  sum(abs(CM.(var1))>=0.1)/size(CM,1)]];

var2='RDI_RScene_field';
dataR = [[sum(abs(SS.(var2))>=0.1)/size(SS,1) sum(abs(SM.(var2))>=0.1)/size(SM,1)],...
    [sum(abs(CS.(var2))>=0.1)/size(CS,1)  sum(abs(CM.(var2))>=0.1)/size(CM,1)]];

var3='RDI_LR_field';
dataC =[[sum(abs(SS.(var3))>=0.1)/size(SS,1) sum(abs(SM.(var3))>=0.1)/size(SM,1)],...
    [sum(abs(CS.(var3))>=0.1)/size(CS,1)  sum(abs(CM.(var3))>=0.1)/size(CM,1)]];

data_all = [sum((abs(SS.(var1))>=0.1) & ~(abs(SS.(var2))>=0.1) & ~(abs(SS.(var3))>=0.1)),...
    sum(~(abs(SS.(var1))>=0.1) & (abs(SS.(var2))>=0.1) & ~(abs(SS.(var3))>=0.1)),...
    sum(~(abs(SS.(var1))>=0.1) & ~(abs(SS.(var2))>=0.1) & (abs(SS.(var3))>=0.1)),...
    sum((abs(SS.(var1))>=0.1) & (abs(SS.(var2))>=0.1) & ~(abs(SS.(var3))>=0.1)),...
    sum((abs(SS.(var1))>=0.1) & ~(abs(SS.(var2))>=0.1) & (abs(SS.(var3))>=0.1)),...
    sum(~(abs(SS.(var1))>=0.1) & (abs(SS.(var2))>=0.1) & (abs(SS.(var3))>=0.1)),...
    sum((abs(SS.(var1))>=0.1) & (abs(SS.(var2))>=0.1) & (abs(SS.(var3))>=0.1)),...
    sum((abs(SS.(var1))>=0.1) | (abs(SS.(var2))>=0.1) | (abs(SS.(var3))>=0.1)),size(SS,1);...

    sum((abs(SM.(var1))>=0.1) & ~(abs(SM.(var2))>=0.1) & ~(abs(SM.(var3))>=0.1)),...
    sum(~(abs(SM.(var1))>=0.1) & (abs(SM.(var2))>=0.1) & ~(abs(SM.(var3))>=0.1)),...
    sum(~(abs(SM.(var1))>=0.1) & ~(abs(SM.(var2))>=0.1) & (abs(SM.(var3))>=0.1)),...
    sum((abs(SM.(var1))>=0.1) & (abs(SM.(var2))>=0.1) & ~(abs(SM.(var3))>=0.1)),...
    sum((abs(SM.(var1))>=0.1) & ~(abs(SM.(var2))>=0.1) & (abs(SM.(var3))>=0.1)),...
    sum(~(abs(SM.(var1))>=0.1) & (abs(SM.(var2))>=0.1) & (abs(SM.(var3))>=0.1)),...
    sum((abs(SM.(var1))>=0.1) & (abs(SM.(var2))>=0.1) & (abs(SM.(var3))>=0.1)),...
    sum((abs(SM.(var1))>=0.1) | (abs(SM.(var2))>=0.1) | (abs(SM.(var3))>=0.1)), size(SM,1);...

    sum((abs(CS.(var1))>=0.1) & ~(abs(CS.(var2))>=0.1) & ~(abs(CS.(var3))>=0.1)),...
    sum(~(abs(CS.(var1))>=0.1) & (abs(CS.(var2))>=0.1) & ~(abs(CS.(var3))>=0.1)),...
    sum(~(abs(CS.(var1))>=0.1) & ~(abs(CS.(var2))>=0.1) & (abs(CS.(var3))>=0.1)),...
    sum((abs(CS.(var1))>=0.1) & (abs(CS.(var2))>=0.1) & ~(abs(CS.(var3))>=0.1)),...
    sum((abs(CS.(var1))>=0.1) & ~(abs(CS.(var2))>=0.1) & (abs(CS.(var3))>=0.1)),...
    sum(~(abs(CS.(var1))>=0.1) & (abs(CS.(var2))>=0.1) & (abs(CS.(var3))>=0.1)),...
    sum((abs(CS.(var1))>=0.1) & (abs(CS.(var2))>=0.1) & (abs(CS.(var3))>=0.1)),...
    sum((abs(CS.(var1))>=0.1) | (abs(CS.(var2))>=0.1) | (abs(CS.(var3))>=0.1)), size(CS,1);...

    sum((abs(CM.(var1))>=0.1) & ~(abs(CM.(var2))>=0.1) & ~(abs(CM.(var3))>=0.1)),...
    sum(~(abs(CM.(var1))>=0.1) & (abs(CM.(var2))>=0.1) & ~(abs(CM.(var3))>=0.1)),...
    sum(~(abs(CM.(var1))>=0.1) & ~(abs(CM.(var2))>=0.1) & (abs(CM.(var3))>=0.1)),...
    sum((abs(CM.(var1))>=0.1) & (abs(CM.(var2))>=0.1) & ~(abs(CM.(var3))>=0.1)),...
    sum((abs(CM.(var1))>=0.1) & ~(abs(CM.(var2))>=0.1) & (abs(CM.(var3))>=0.1)),...
    sum(~(abs(CM.(var1))>=0.1) & (abs(CM.(var2))>=0.1) & (abs(CM.(var3))>=0.1)),...
    sum((abs(CM.(var1))>=0.1) & (abs(CM.(var2))>=0.1) & (abs(CM.(var3))>=0.1)),...
    sum((abs(CM.(var1))>=0.1) | (abs(CM.(var2))>=0.1) | (abs(CM.(var3))>=0.1)), size(CM,1)];


figure;
subplot(2,2,1); hold on
dat=[];
for i=1:4
dat(i,1) = (data_all(i,1)+data_all(i,2)+data_all(i,4))/data_all(i,9);
dat(i,2) = (data_all(i,3))/data_all(i,9);
dat(i,3) = (data_all(i,5)+data_all(i,6)+data_all(i,7))/data_all(i,9);
end
bar([1 1 1 1]); bar(dat,'stacked')
legend({'all','Scene','Choice','S & C'})

subplot(2,2,2); hold on
dat=[];
for i=1:4
dat(i,1) = (data_all(i,1))/data_all(i,9);
dat(i,2) = (data_all(i,4))/data_all(i,9);
dat(i,3) = (data_all(i,5))/data_all(i,9);
dat(i,4) = (data_all(i,7))/data_all(i,9);
end
bar([1 1 1 1]); bar(dat,'stacked')
legend({'all','Left','L&R','L&C', 'L&R&C'})

subplot(2,2,3); hold on
dat=[];
for i=1:4
dat(i,1) = (data_all(i,2))/data_all(i,9);
dat(i,2) = (data_all(i,5))/data_all(i,9);
dat(i,3) = (data_all(i,6))/data_all(i,9);
dat(i,4) = (data_all(i,7))/data_all(i,9);
end
bar([1 1 1 1]); bar(dat,'stacked')
legend({'all','right','L&R','R&C', 'L&R&C'})

subplot(2,2,4); hold on
dat=[];
for i=1:4
dat(i,1) = (data_all(i,3))/data_all(i,9);
dat(i,2) = (data_all(i,4))/data_all(i,9);
dat(i,3) = (data_all(i,6))/data_all(i,9);
dat(i,4) = (data_all(i,7))/data_all(i,9);
end
bar([1 1 1 1]); bar(dat,'stacked')
legend({'all','Choice','L&C','R&C', 'L&R&C'})
%%

    X.unit(x) = X{x}
[C,ia] = unique(X.rat + 0.001*X.AvgFR);
B = X(ia,:)

X = UnF_SUB(max(abs([UnF_SUB.RDI_LScene,UnF_SUB.RDI_RScene,UnF_SUB.RDI_LR]),[],2)>=0.1,:)
[C,ia] = unique(X.rat + 0.01*X.session+0.0001*X.TT+0.000001*X.ID{end-1:end});
C = X(ia,:)
%% 2 items / 4 items box plot
var = 'onMazeAvgFR';
var_label = 'onMazeAvgFR';
x=Un_SUB.(var); y=Un_CA1.(var);
x1 = x(Un_SUB.NumField==1); x2 = x(Un_SUB.NumField>1); 
y1 = y(Un_CA1.NumField==1); y2 = y(Un_CA1.NumField>1); 
data1 = [[x; nan(length(y)-length(x),1)],y];
data2 = [[x1;nan(length(y1)-length(x1),1)],[x2;nan(length(y1)-length(x2),1)],y1,[y2;nan(length(y1)-length(y2),1)]];
% [~,p1,~,stats1] = ttest2(x,y);
% [~,p2,~,stats2] = ttest2(x1,x2);
% [~,p3,~,stats3] = ttest2(y1,y2);
% [~,p4,~,stats4] = ttest2(x1,y1);
% [~,p5,~,stats5] = ttest2(x2,y2);

[p1,h,stats1] = ranksum(x,y);
[~,p2,~,stats2] = ttest2(x1,x2)
[~,p3,~,stats3] = ranksum(y1,y2)
[~,p4,~,stats4] = ttest2(x1,y1)
[~,p5,~,stats5] = ttest2(x2,y2)


figure;
boxplot(data1,'symbol','o')
title(['t=' jjnum2str(stats1.zval,3) ',p=' jjnum2str(p1,3)])
ylabel(var_label)
xticklabels({'SUB','CA1'})
saveas(gca,[ROOT.Save '\Manuscript\fig2\' var '_region.png'])
saveas(gca,[ROOT.Save '\Manuscript\fig2\' var '_region.svg'])


% figure;
% boxplot(data2,'symbol','o')
% title(['t=' jjnum2str(stats2.zval,3) ',' jjnum2str(stats3.tstat,3) ',' jjnum2str(stats4.tstat,3) ',' jjnum2str(stats5.tstat,3) ...
%     ',p=' jjnum2str(p2,3) ',' jjnum2str(p3,3) ',' jjnum2str(p4,3) ',' jjnum2str(p5,3)])
% ylabel(var_label)
% xticklabels({'SUB_SF','SUB_MF','CA1_SF','CA1_MF'})
% saveas(gca,[ROOT.Save '\Manuscript\fig2\' var '_fields.png'])
% saveas(gca,[ROOT.Save '\Manuscript\fig2\' var '_fields.svg'])


%% add nFields
% 
% Un = Un_CA1; UnF=UnF_CA1;
% for uid = 1:size(UnF,1)
%     UID = UnF.ID{uid};
%     UID_o = strncmp(UID,Un.ID,12);
%     UnF.NumField(uid) = Un.NumField(UID_o);
% 
%        Dat = load([ROOT.Save '\units_mat\U1\' UID '.mat'],'RDIs_field','thisFieldMap1D');
% 
%        UnF.onMazeAvgFR_field(uid) = Dat.thisFieldMap1D.onmazeAvgFR1D(1);
%        UnF.SI_field(uid) = Dat.thisFieldMap1D.SpaInfoScore1D(1);
% 
% end
% writetable(UnF,[ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx'],'writemode','replacefile');
% 
rdis = [(UnF_SUB.RDI_LScene),(UnF_SUB.RDI_RScene),(UnF_SUB.RDI_LR)];
[m,t] = nanmax([abs(UnF_SUB.RDI_LScene),abs(UnF_SUB.RDI_RScene),abs(UnF_SUB.RDI_LR)],[],2);
for u=1:size(UnF_SUB,1)
UnF_SUB.RDI_Max(u) = rdis(u,t(u));
end

rdis = [(UnF_CA1.RDI_LScene),(UnF_CA1.RDI_RScene),(UnF_CA1.RDI_LR)];
[m,t] = nanmax([abs(UnF_CA1.RDI_LScene),abs(UnF_CA1.RDI_RScene),abs(UnF_CA1.RDI_LR)],[],2);
for u=1:size(UnF_CA1,1)
UnF_CA1.RDI_Max(u) = rdis(u,t(u));
end


for u=1:size(Un_SUB,1)
    idx = find(strncmp(Un_SUB.ID{u},UnF_SUB.ID,12));
Un_SUB.RDI_Max(u) = nanmax(abs(UnF_SUB.RDI_Max(idx)));
end

for u=1:size(Un_CA1,1)
    idx = find(strncmp(Un_CA1.ID{u},UnF_CA1.ID,12));
Un_CA1.RDI_Max(u) = nanmax(abs(UnF_CA1.RDI_Max(idx)));
end

%% 2 items / 4 items bar plot _ TP fields
var = 'onMazeAvgFR_field';
var_label = 'mean firing rate (Hz)';
x=UnF_SUB.(var); y=UnF_CA1.(var);
x1 = x(UnF_SUB.NumField==1); x2 = x(UnF_SUB.NumField>1); 
y1 = y(UnF_CA1.NumField==1); y2 = y(UnF_CA1.NumField>1); 
m1 = max([length(x),length(y)]); m2 = max([length(x1),length(x2),length(y1),length(y2)]);
data1 = [[x; nan(m1-length(x),1)],[y; nan(m1-length(y),1)]];
data2 = [[x1;nan(m2-length(x1),1)],[x2;nan(m2-length(x2),1)],...
    [y1; nan(m2-length(y1),1)],[y2;nan(m2-length(y2),1)]];
[~,p1,~,stats1] = ttest2(x,y);
[~,p2,~,stats2] = ttest2(x1,x2);
[~,p3,~,stats3] = ttest2(y1,y2);
[~,p4,~,stats4] = ttest2(x1,y1);
[~,p5,~,stats5] = ttest2(x2,y2);

data = [x1;x2;y1;y2];
g1 = [ones(size([x1;x2],1),1);2*ones(size([y1;y2],1),1)];
g2 = [ones(size([x1],1),1);2*ones(size([x2],1),1);ones(size([y1],1),1);2*ones(size([y2],1),1)];
anovan(data,{g1 g2},'model','interaction','varnames',{'field','region'})



figure; hold on
dat = nanmean(data1);
err = nanstd(data1) ./ [sqrt(length(x)) sqrt(length(y))];
bar(dat)
errorbar(dat,err)
title(['t=' jjnum2str(stats1.tstat,3) ',p=' jjnum2str(p1,3)])
ylabel(var_label)
xticklabels({'SUB','CA1'})
% saveas(gca,[ROOT.Mother '\Manuscript Figures\fig3\(bar)' var '_region.png'])
% saveas(gca,[ROOT.Mother '\Manuscript Figures\fig3\(bar)' var '_region.svg'])


figure; hold on
dat = nanmean(data2);
err = nanstd(data2) ./ [sqrt(length(x1)) sqrt(length(x2)) sqrt(length(y1)) sqrt(length(y2))];
bar(dat)
errorbar(dat,err)

title(['t=' jjnum2str(stats2.tstat,3) ',' jjnum2str(stats3.tstat,3) ',' jjnum2str(stats4.tstat,3) ',' jjnum2str(stats5.tstat,3) ...
    ',p=' jjnum2str(p2,3) ',' jjnum2str(p3,3) ',' jjnum2str(p4,3) ',' jjnum2str(p5,3)])
ylabel(var_label)
xticklabels({'SUB_SF','SUB_MF','CA1_SF','CA1_MF'})
% saveas(gca,[ROOT.Mother '\Manuscript Figures\fig3\(bar)' var '_fields.png'])
% saveas(gca,[ROOT.Mother '\Manuscript Figures\fig3\(bar)' var '_fields.svg'])
%% 2 items / 4 items box plot _ TP fields
var = 'RDI_Max';
var_label = 'Selectivity Index (Max)';
x=UnF_SUB.(var); y=UnF_CA1.(var);
x1 = x(UnF_SUB.NumField==1); x2 = x(UnF_SUB.NumField>1); 
y1 = y(UnF_CA1.NumField==1); y2 = y(UnF_CA1.NumField>1); 
data1 = [[x],[y; nan(length(x)-length(y),1)]];
data2 = [[x1;nan(length(x2)-length(x1),1)],[x2;nan(length(x2)-length(x2),1)],...
    [y1; nan(length(x2)-length(y1),1)],[y2;nan(length(x2)-length(y2),1)]];
[~,p1,~,stats1] = ttest2(x,y);
[~,p2,~,stats2] = ttest2(x1,x2);
[~,p3,~,stats3] = ttest2(y1,y2);
[~,p4,~,stats4] = ttest2(x1,y1);
[~,p5,~,stats5] = ttest2(x2,y2);

[~,k1] = kstest2(x,y);
[~,k2] = kstest2(x1,x2);
[~,k3] = kstest2(y1,y2);
[~,k4] = kstest2(x1,y1);
[~,k5] = kstest2(x2,y2);

figure;
boxplot(data1,'symbol','.')
title(['t=' jjnum2str(stats1.tstat,3) ',p=' jjnum2str(p1,3) ',k=' jjnum2str(k1,3)])
ylabel(var_label)
xticklabels({'SUB','CA1'})
saveas(gca,[ROOT.Mother '\Manuscript Figures\fig3\' var '_region_TP.png'])
saveas(gca,[ROOT.Mother '\Manuscript Figures\fig3\' var '_region_TP.svg'])


figure;
boxplot(data2,'symbol','.')
title(['t=' jjnum2str(stats2.tstat,3) ',' jjnum2str(stats3.tstat,3) ',' jjnum2str(stats4.tstat,3) ',' jjnum2str(stats5.tstat,3) ...
    ',p=' jjnum2str(p2,3) ',' jjnum2str(p3,3) ',' jjnum2str(p4,3) ',' jjnum2str(p5,3)...
     ',k=' jjnum2str(k2,3) ',' jjnum2str(k3,3) ',' jjnum2str(k4,3) ',' jjnum2str(k5,3)])
ylabel(var_label)
xticklabels({'SUB_SF','SUB_MF','CA1_SF','CA1_MF'})
saveas(gca,[ROOT.Mother '\Manuscript Figures\fig3\' var '_fields_TP.png'])
saveas(gca,[ROOT.Mother '\Manuscript Figures\fig3\' var '_fields_TP.svg'])
%% numfields_pie
var = 'NumField';
var_label = 'Number of fields';
x=Un_SUB.(var); y=Un_CA1.(var);

x1=[]; y1=[];
for i=1:5
x1(i) = sum(x==i); y1(i) = sum(y==i);
end

figure;
hold on
subplot(1,2,1)
pie(flip(x1))
title(['SUB ' num2str(x1)])
subplot(1,2,2)
pie(flip(y1))
title(['CA1 ' num2str(y1)])

%% numfields_pie
var = 'Selectivity_LScene';
var_label = 'Unit type (LScene)';
x=Un_SUB.(var); y=Un_CA1.(var);

x1=[]; y1=[];
for i=0:8
x1(i+1) = sum(x==i); y1(i+1) = sum(y==i);
end

figure;
hold on
subplot(1,2,1)
pie(flip(x1))
title('SUB')
subplot(1,2,2)
pie(flip(y1))
title('CA1')

%% cdf_RDI
figure
varList = {'LScene','RScene','LR'};
for v=1:3
subplot(2,3,v); hold on
c0 = cdfplot(abs(UnF_SUB.(['RDI_' varList{v}])));
c0.Color = CList(1,:);
c1 = cdfplot(abs(UnF_CA1.(['RDI_' varList{v}])));
c1.Color = CList(2,:);
% [h,p] = kstest2((R0.(['m' var])),(R1.(['m' var])))
xlim([0 1.2])


subplot(2,3,v+3); hold on
x0 = abs(UnF_SUB.(['RDI_' varList{v}]));
x1 = abs(UnF_CA1.(['RDI_' varList{v}]));
dat = [nanmean(x0) nanmean(x1)];
err = [nanstd(x0)/sqrt(sum(~isnan(x0))) nanstd(x1)/sqrt(sum(~isnan(x1)))];
bar(dat)
errorbar(dat,err)
% yticks([0:0.05:0.2])
[p,h,stats] = ranksum(x0,x1)
title(['p=' jjnum2str(p,3) ', t=' jjnum2str(stats.zval,3)])


end

%% hist_RDI
figure
varList = {'LScene','RScene','LR'};
for v=1:3
subplot(2,3,v); hold on
c0 = histogram((UnF_SUB.(['RDI_' varList{v}])),'binwidth',.1,'Normalization','probability');
c0.FaceColor = CList(1,:);
xlim([-1.2 1.2]);ylim([0 .2])

subplot(2,3,v+3); hold on
c1 = histogram((UnF_CA1.(['RDI_' varList{v}])),'binwidth',.1,'Normalization','probability');
c1.FaceColor = CList(2,:);
xlim([-1.2 1.2]);ylim([0 .2])

[h,p,r] = kstest2((UnF_CA1.(['RDI_' varList{v}])),(UnF_SUB.(['RDI_' varList{v}])))


title(['p = ' num2str(p)])




end