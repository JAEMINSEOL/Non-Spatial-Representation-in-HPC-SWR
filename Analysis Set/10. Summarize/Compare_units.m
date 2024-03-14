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
UnitsTable.SUB_FR = readtable([ROOT.Save '\UnitsTable_SUB_RDI_FR.xlsx']);
UnitsTable_field.SUB = readtable([ROOT.Save '\UnitsTable_SUB_field_forAnalysis.xlsx']);
% UnitPair.SUB = readtable([ROOT.Save '\UnitPair_SUB.xlsx']);
% UnitPair_field.SUB = readtable([ROOT.Save '\UnitPair_SUB_field.xlsx']);

RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx']);
ReactTable.CA1 = readtable([ROOT.Save '\ReactTable_CA1_CA1.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_forAnalysis_TP.xlsx']);
UnitsTable.CA1_FR = readtable([ROOT.Save '\UnitsTable_CA1_RDI_FR.xlsx']);
UnitsTable_field.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx']);
% UnitPair.CA1= readtable([ROOT.Save '\UnitPair_CA1.xlsx']);
% UnitPair_field.CA1 = readtable([ROOT.Save '\UnitPair_CA1_field.xlsx']);

Un_SUB = UnitsTable.SUB;
Un_CA1= UnitsTable.CA1;
UnF_SUB = UnitsTable_field.SUB;
UnF_CA1= UnitsTable_field.CA1;

Un_SUB_FR = UnitsTable.SUB_FR;
Un_CA1_FR= UnitsTable.CA1_FR;

Un_SUB(Un_SUB.rat==232 & Un_SUB.session==4,:)=[];
UnF_SUB(UnF_SUB.rat==232 & UnF_SUB.session==4,:)=[];
Un_CA1(Un_CA1.rat==232 & Un_CA1.session==4,:)=[];
UnF_CA1(UnF_CA1.rat==232 & UnF_CA1.session==4,:)=[];
Un_SUB_FR(Un_SUB_FR.rat==232 & Un_SUB_FR.session==4,:)=[];
Un_CA1_FR(Un_CA1_FR.rat==232 & Un_CA1_FR.session==4,:)=[];


for u=1:size(Un_CA1)
    tFields = UnF_CA1(find(strncmp(Un_CA1.ID{u},UnF_CA1.ID,12)),:);
    Un_CA1.nFields(u)=size(tFields,1);
    if size(tFields,1)<1
        continue;
    end
    [~,t] = max(abs(tFields.RDI_LScene)); Un_CA1.RDI_LScene(u)=tFields.RDI_LScene(t);
    [~,t] = max(abs(tFields.RDI_RScene)); Un_CA1.RDI_RScene(u)=tFields.RDI_RScene(t);
    [~,t] = max(abs(tFields.RDI_LR)); Un_CA1.RDI_LR(u)=tFields.RDI_LR(t);
end

for u=1:size(Un_SUB)
    tFields = UnF_SUB(find(strncmp(Un_SUB.ID{u},UnF_SUB.ID,12)),:);
    Un_SUB.nFields(u)=size(tFields,1);
    if size(tFields,1)<1
        continue;
    end
    [~,t] = max(abs(tFields.RDI_LScene)); Un_SUB.RDI_LScene(u)=tFields.RDI_LScene(t);
    [~,t] = max(abs(tFields.RDI_RScene)); Un_SUB.RDI_RScene(u)=tFields.RDI_RScene(t);
    [~,t] = max(abs(tFields.RDI_LR)); Un_SUB.RDI_LR(u)=tFields.RDI_LR(t);
end

%% histogram_RDI_FR vs. TP
sz=0.1;
figure;
subplot(3,2,1); hold on; ylim([0 .2]); xlim([-1.2 1.2]); title('SUB')
histogram(Un_SUB.RDI_LScene,'BinWidth',sz,'Normalization','probability');
histogram(Un_SUB_FR.RDI_LScene,'BinWidth',sz,'Normalization','probability');
subplot(3,2,2); hold on; ylim([0 .2]); xlim([-1.2 1.2]); title('CA1')
histogram(Un_CA1.RDI_LScene,'BinWidth',sz,'Normalization','probability');
histogram(Un_CA1_FR.RDI_LScene,'BinWidth',sz,'Normalization','probability');

subplot(3,2,3); hold on; ylim([0 .2]); xlim([-1.2 1.2])
histogram(Un_SUB.RDI_RScene,'BinWidth',sz,'Normalization','probability');
histogram(Un_SUB_FR.RDI_RScene,'BinWidth',sz,'Normalization','probability');
subplot(3,2,4); hold on; ylim([0 .2])
histogram(Un_CA1.RDI_RScene,'BinWidth',sz,'Normalization','probability');
histogram(Un_CA1_FR.RDI_RScene,'BinWidth',sz,'Normalization','probability');

subplot(3,2,5); hold on; ylim([0 .2]); xlim([-1.2 1.2])
histogram(Un_SUB.RDI_LR,'BinWidth',sz,'Normalization','probability');
histogram(Un_SUB_FR.RDI_LR,'BinWidth',sz,'Normalization','probability');
subplot(3,2,6); hold on; ylim([0 .2]); xlim([-1.2 1.2])
histogram(Un_CA1.RDI_LR,'BinWidth',sz,'Normalization','probability');
histogram(Un_CA1_FR.RDI_LR,'BinWidth',sz,'Normalization','probability');

%% histogram_RDI_MF vs. SF
figure;
 subplot(3,2,1);hold on; ylim([0 .2]); xlim([-1.2 1.2]); title(['n=' num2str(sum(Un_SUB.nFields>1)) ',' num2str(sum(Un_SUB.nFields==1))])
histogram(Un_SUB.RDI_LScene(Un_SUB.nFields>1),'BinWidth',sz,'Normalization','probability');
histogram(Un_SUB.RDI_LScene(Un_SUB.nFields==1),'BinWidth',sz,'Normalization','probability');
subplot(3,2,2); hold on; ylim([0 .2]); xlim([-1.2 1.2]); title(['n=' num2str(sum(Un_CA1.nFields>1)) ',' num2str(sum(Un_CA1.nFields==1))])
histogram(Un_CA1.RDI_LScene(Un_CA1.nFields>1),'BinWidth',sz,'Normalization','probability');
histogram(Un_CA1.RDI_LScene(Un_CA1.nFields==1),'BinWidth',sz,'Normalization','probability');

 subplot(3,2,3);hold on; ylim([0 .2]); xlim([-1.2 1.2])
histogram(Un_SUB.RDI_RScene(Un_SUB.nFields>1),'BinWidth',sz,'Normalization','probability');
histogram(Un_SUB.RDI_RScene(Un_SUB.nFields==1),'BinWidth',sz,'Normalization','probability');
subplot(3,2,4); hold on; ylim([0 .2]); xlim([-1.2 1.2])
histogram(Un_CA1.RDI_RScene(Un_CA1.nFields>1),'BinWidth',sz,'Normalization','probability');
histogram(Un_CA1.RDI_RScene(Un_CA1.nFields==1),'BinWidth',sz,'Normalization','probability');

 subplot(3,2,5);hold on; ylim([0 .2]); xlim([-1.2 1.2])
histogram(Un_SUB.RDI_LR(Un_SUB.nFields>1),'BinWidth',sz,'Normalization','probability');
histogram(Un_SUB.RDI_LR(Un_SUB.nFields==1),'BinWidth',sz,'Normalization','probability');
subplot(3,2,6); hold on; ylim([0 .2]); xlim([-1.2 1.2])
histogram(Un_CA1.RDI_LR(Un_CA1.nFields>1),'BinWidth',sz,'Normalization','probability');
histogram(Un_CA1.RDI_LR(Un_CA1.nFields==1),'BinWidth',sz,'Normalization','probability');

%% bar_RDI_MF vs SF
vsList={'RDI_LScene','RDI_RScene','RDI_LR'};
U0 = Un_SUB;
U1 = Un_CA1;
US0 = U0(U0.NumField==1,:);
US1 = U0(U0.NumField>1,:);
UC0 = U1(U1.NumField==1,:);
UC1 = U1(U1.NumField>1,:);

xS0L = abs(US0.RDI_LScene); xS1L = abs(US1.RDI_LScene);
xC0L = abs(UC0.RDI_LScene); xC1L = abs(UC1.RDI_LScene);
xS0R = abs(US0.RDI_RScene); xS1R = abs(US1.RDI_RScene);
xC0R = abs(UC0.RDI_RScene); xC1R = abs(UC1.RDI_RScene);
xS0C = abs(US0.RDI_LR); xS1C = abs(US1.RDI_LR);
xC0C = abs(UC0.RDI_LR); xC1C = abs(UC1.RDI_LR);

figure; 
subplot(1,2,1); hold on
dat = [nanmean(xS0L) nanmean(xS1L) nanmean(xS0R) nanmean(xS1R) nanmean(xS0C) nanmean(xS1C)];
err = [nanstd(xS0L)/sqrt(sum(~isnan(xS0L))) nanstd(xS1L)/sqrt(sum(~isnan(xS1L)))...
    nanstd(xS0R)/sqrt(sum(~isnan(xS0R))) nanstd(xS1R)/sqrt(sum(~isnan(xS1R)))...
    nanstd(xS0C)/sqrt(sum(~isnan(xS0C))) nanstd(xS1C)/sqrt(sum(~isnan(xS1C)))];
bar(dat)
errorbar(dat,err)
ylim([0 .5])
title('RDI SF vs. MF, SUB')
xticks([1.5:2:5.5]); xticklabels({'Left', 'Right', 'Choice'})

subplot(1,2,2); hold on
dat = [nanmean(xC0L) nanmean(xC1L) nanmean(xC0R) nanmean(xC1R) nanmean(xC0C) nanmean(xC1C)];
err = [nanstd(xC0L)/sqrt(sum(~isnan(xC0L))) nanstd(xC1L)/sqrt(sum(~isnan(xC1L)))...
    nanstd(xC0R)/sqrt(sum(~isnan(xC0R))) nanstd(xC1R)/sqrt(sum(~isnan(xC1R)))...
    nanstd(xC0C)/sqrt(sum(~isnan(xC0C))) nanstd(xC1C)/sqrt(sum(~isnan(xC1C)))];
bar(dat)
errorbar(dat,err)
ylim([0 .5])
title('RDI SF vs. MF, CA1')
xticks([1.5:2:5.5]); xticklabels({'Left', 'Right', 'Choice'})
% % yticks([0:al:0.2])

tab_RDI=table;
[p,h,stats] = ranksum(xS0L,xS1L);
tab_RDI.SUB_Left(1)=p; tab_RDI.SUB_Left(2)=stats.zval;
[p,h,stats] = ranksum(xS0R,xS1R);
tab_RDI.SUB_Right(1)=p; tab_RDI.SUB_Right(2)=stats.zval;
[p,h,stats] = ranksum(xS0C,xS1C);
tab_RDI.SUB_Choice(1)=p; tab_RDI.SUB_Choice(2)=stats.zval;
[p,h,stats] = ranksum(xC0L,xC1L);
tab_RDI.CA1_Left(1)=p; tab_RDI.CA1_Left(2)=stats.zval;
[p,h,stats] = ranksum(xC0R,xC1R);
tab_RDI.CA1_Right(1)=p; tab_RDI.CA1_Right(2)=stats.zval;
[p,h,stats] = ranksum(xC0C,xC1C);
tab_RDI.CA1_Choice(1)=p; tab_RDI.CA1_Choice(2)=stats.zval;


% ANOVA
data = [xS0L;xS0R;xS0C;xC0L;xC0R;xC0C;xS1L;xS1R;xS1C;xC1L;xC1R;xC1C];
g1 = [repmat('SUB',[size([xS0L;xS0R;xS0C],1) 1]); repmat('CA1',[size([xC0L;xC0R;xC0C],1) 1]);...
    repmat('SUB',[size([xS1L;xS1R;xS1C],1) 1]); repmat('CA1',[size([xC1L;xC1R;xC1C],1) 1])];
g2 = [repmat('L',[size([xS0L],1) 1]); repmat('R',[size([xS0R],1) 1]); repmat('C',[size([xS0C],1) 1]);...
    repmat('L',[size([xC0L],1) 1]); repmat('R',[size([xC0R],1) 1]); repmat('C',[size([xC0C],1) 1]);...
    repmat('L',[size([xS1L],1) 1]); repmat('R',[size([xS1R],1) 1]); repmat('C',[size([xS1C],1) 1]);...
    repmat('L',[size([xC1L],1) 1]); repmat('R',[size([xC1R],1) 1]); repmat('C',[size([xC1C],1) 1]);];
g3 = [repmat('SF',[size([xS0L],1) 1]); repmat('SF',[size([xS0R],1) 1]); repmat('SF',[size([xS0C],1) 1]);...
    repmat('SF',[size([xC0L],1) 1]); repmat('SF',[size([xC0R],1) 1]); repmat('SF',[size([xC0C],1) 1]);...
    repmat('MF',[size([xS1L],1) 1]); repmat('MF',[size([xS1R],1) 1]); repmat('MF',[size([xS1C],1) 1]);...
    repmat('MF',[size([xC1L],1) 1]); repmat('MF',[size([xC1R],1) 1]); repmat('MF',[size([xC1C],1) 1]);];

dat_xls=table;

dat_xls.Region=g1;
dat_xls.Item=g2;
dat_xls.Field=g3;
dat_xls.RPR=data;

writetable(dat_xls,[ROOT.Processed '\SI_SFMF.xlsx'],'writemode','replacefile')

[p,tbl,stats] = anovan(data, {g1 g3},'model', 'interaction');
%% bar_RDI_FR vs. TP
vsList={'RDI_LScene','RDI_RScene','RDI_LR'};

US0 = Un_SUB_FR;
US1 = Un_SUB;
UC0 = Un_CA1_FR;
UC1 = Un_CA1;

xS0L = abs(US0.RDI_LScene); xS1L = abs(US1.RDI_LScene);
xC0L = abs(UC0.RDI_LScene); xC1L = abs(UC1.RDI_LScene);
xS0R = abs(US0.RDI_RScene); xS1R = abs(US1.RDI_RScene);
xC0R = abs(UC0.RDI_RScene); xC1R = abs(UC1.RDI_RScene);
xS0C = abs(US0.RDI_LR); xS1C = abs(US1.RDI_LR);
xC0C = abs(UC0.RDI_LR); xC1C = abs(UC1.RDI_LR);

figure; 
subplot(1,2,1); hold on
dat = [nanmean(xS0L) nanmean(xS1L) nanmean(xS0R) nanmean(xS1R) nanmean(xS0C) nanmean(xS1C)];
err = [nanstd(xS0L)/sqrt(sum(~isnan(xS0L))) nanstd(xS1L)/sqrt(sum(~isnan(xS1L)))...
    nanstd(xS0R)/sqrt(sum(~isnan(xS0R))) nanstd(xS1R)/sqrt(sum(~isnan(xS1R)))...
    nanstd(xS0C)/sqrt(sum(~isnan(xS0C))) nanstd(xS1C)/sqrt(sum(~isnan(xS1C)))];
bar(dat)
errorbar(dat,err)
ylim([0 .5])
title('RDI FR vs. TP, SUB')
xticks([1.5:2:5.5]); xticklabels({'Left', 'Right', 'Choice'})

subplot(1,2,2); hold on
dat = [nanmean(xC0L) nanmean(xC1L) nanmean(xC0R) nanmean(xC1R) nanmean(xC0C) nanmean(xC1C)];
err = [nanstd(xC0L)/sqrt(sum(~isnan(xC0L))) nanstd(xC1L)/sqrt(sum(~isnan(xC1L)))...
    nanstd(xC0R)/sqrt(sum(~isnan(xC0R))) nanstd(xC1R)/sqrt(sum(~isnan(xC1R)))...
    nanstd(xC0C)/sqrt(sum(~isnan(xC0C))) nanstd(xC1C)/sqrt(sum(~isnan(xC1C)))];
bar(dat)
errorbar(dat,err)
ylim([0 .5])
title('RDI FR vs. TP, CA1')
xticks([1.5:2:5.5]); xticklabels({'Left', 'Right', 'Choice'})
% % yticks([0:al:0.2])

tab_RDI=table;
[p,h,stats] = ranksum(xS0L,xS1L);
tab_RDI.SUB_Left(1)=p; tab_RDI.SUB_Left(2)=stats.zval;
[p,h,stats] = ranksum(xS0R,xS1R);
tab_RDI.SUB_Right(1)=p; tab_RDI.SUB_Right(2)=stats.zval;
[p,h,stats] = ranksum(xS0C,xS1C);
tab_RDI.SUB_Choice(1)=p; tab_RDI.SUB_Choice(2)=stats.zval;
[p,h,stats] = ranksum(xC0L,xC1L);
tab_RDI.CA1_Left(1)=p; tab_RDI.CA1_Left(2)=stats.zval;
[p,h,stats] = ranksum(xC0R,xC1R);
tab_RDI.CA1_Right(1)=p; tab_RDI.CA1_Right(2)=stats.zval;
[p,h,stats] = ranksum(xC0C,xC1C);
tab_RDI.CA1_Choice(1)=p; tab_RDI.CA1_Choice(2)=stats.zval;

% ANOVA
data = [xS0L;xS0R;xS0C;xC0L;xC0R;xC0C;xS1L;xS1R;xS1C;xC1L;xC1R;xC1C];
g1 = [repmat('SUB',[size([xS0L;xS0R;xS0C],1) 1]); repmat('CA1',[size([xC0L;xC0R;xC0C],1) 1]);...
    repmat('SUB',[size([xS1L;xS1R;xS1C],1) 1]); repmat('CA1',[size([xC1L;xC1R;xC1C],1) 1])];
g2 = [repmat('L',[size([xS0L],1) 1]); repmat('R',[size([xS0R],1) 1]); repmat('C',[size([xS0C],1) 1]);...
    repmat('L',[size([xC0L],1) 1]); repmat('R',[size([xC0R],1) 1]); repmat('C',[size([xC0C],1) 1]);...
    repmat('L',[size([xS1L],1) 1]); repmat('R',[size([xS1R],1) 1]); repmat('C',[size([xS1C],1) 1]);...
    repmat('L',[size([xC1L],1) 1]); repmat('R',[size([xC1R],1) 1]); repmat('C',[size([xC1C],1) 1]);];
g3 = [repmat('FR',[size([xS0L],1) 1]); repmat('FR',[size([xS0R],1) 1]); repmat('FR',[size([xS0C],1) 1]);...
    repmat('FR',[size([xC0L],1) 1]); repmat('FR',[size([xC0R],1) 1]); repmat('FR',[size([xC0C],1) 1]);...
    repmat('TP',[size([xS1L],1) 1]); repmat('TP',[size([xS1R],1) 1]); repmat('TP',[size([xS1C],1) 1]);...
    repmat('TP',[size([xC1L],1) 1]); repmat('TP',[size([xC1R],1) 1]); repmat('TP',[size([xC1C],1) 1]);];

dat_xls=table;

dat_xls.Region=g1;
dat_xls.Item=g2;
dat_xls.Field=g3;
dat_xls.RPR=data;

writetable(dat_xls,[ROOT.Processed '\SI_FRTP.xlsx'])

[p,tbl,stats] = anovan(data, {g1 g3},'model', 'interaction');



%% pie_Unit field heterogeniety
% var = DecodingP_all;
U0 = Un_SUB(isnan(Un_SUB.FieldHet) & Un_SUB.NumField>1,:);
U1 = Un_CA1(~isnan(Un_CA1.FieldHet) & Un_CA1.NumField>1,:);

U0x = U0(isnan(U0.FieldHet),:);   U00 = U0(U0.FieldHet==0,:); U01 = U0(U0.FieldHet==1,:); U02 = U0(U0.FieldHet==2,:); U03 = U0(U0.FieldHet==3,:); 
U1x = U1(isnan(U1.FieldHet),:);   U10 = U1(U1.FieldHet==0,:); U11 = U1(U1.FieldHet==1,:); U12 = U1(U1.FieldHet==2,:); U13 = U1(U1.FieldHet==3,:); 

figure;
subplot(1,2,1)
pie([size(U0x,1), size(U00,1),size(U01,1),size(U02,1),size(U03,1)],'%.3f%%')
title(['SUB - ' num2str(size(U0x,1)) ', ' num2str(size(U00,1)) ', ' num2str(size(U01,1)) ', ' num2str(size(U02,1)) ', ' num2str(size(U03,1))])
subplot(1,2,2)
pie([size(U1x,1), size(U10,1),size(U11,1),size(U12,1),size(U13,1)],'%.3f%%')
title(['CA1 - ' num2str(size(U1x,1)) ', ' num2str(size(U10,1)) ', ' num2str(size(U11,1)) ', ' num2str(size(U12,1)) ', ' num2str(size(U13,1))])

legend({'all fields are Non-sig', 'homo' , 'hetero-scenes','hetero-scene/choice','hetero-all 3'},'location','eastoutside')
%% bar_RPR_het

U00 = U0(U0.FieldHet<=1,:);
U01 = U0(U0.FieldHet>1,:);
U10 = U1(U1.FieldHet<=1,:);
U11 = U1(U1.FieldHet>1,:);

figure; 
subplot(1,2,1); hold on
x00 = U00.RipPartRate_all; x01 = U01.RipPartRate_all;
x10 = U10.RipPartRate_all; x11 = U11.RipPartRate_all;

dat = [nanmean(x00) nanmean(x01) nanmean(x10) nanmean(x11)];
err = [nanstd(x00)/sqrt(sum(~isnan(x00))) nanstd(x01)/sqrt(sum(~isnan(x01))) nanstd(x10)/sqrt(sum(~isnan(x10))) nanstd(x11)/sqrt(sum(~isnan(x11)))];
bar(dat)
errorbar(dat,err)
ylim([0 .6])
title('Ripple Participation Rate, all SWRs')
xticks([1:4]); xticklabels({'SUB-Hom','SUB-Het','CA1-Hom','CA1-Het'})

subplot(1,2,2); hold on
x00 = U00.RipPartRate_NonSp; x01 = U01.RipPartRate_NonSp;
x10 = U10.RipPartRate_NonSp; x11 = U11.RipPartRate_NonSp;

dat = [nanmean(x00) nanmean(x01) nanmean(x10) nanmean(x11)];
err = [nanstd(x00)/sqrt(sum(~isnan(x00))) nanstd(x01)/sqrt(sum(~isnan(x01))) nanstd(x10)/sqrt(sum(~isnan(x10))) nanstd(x11)/sqrt(sum(~isnan(x11)))];
bar(dat)
errorbar(dat,err)
ylim([0 .6])
title('Ripple Participation Rate, task-related SWRs')
xticks([1:4]); xticklabels({'SUB-Hom','SUB-Het','CA1-Hom','CA1-Het'})
% % yticks([0:al:0.2])
[p,h,stats] = ranksum(x01,x11)
title(['p=' jjnum2str(p,3) ', t=' jjnum2str(stats.zval,3)])

%%
data = [x00;x01;x10;x11];
g1 = [repmat('SUB',[size([x00;x01],1) 1]); repmat('CA1',[size([x10;x11],1) 1])]
g2 = [repmat('Het',[size([x00],1) 1]); repmat('Hom',[size([x01],1) 1]);repmat('Het',[size([x10],1) 1]); repmat('Hom',[size([x11],1) 1])];

[p,tbl,stats] = anovan(data, {g1,g2}, 'model', 'interaction');

%% bar_RPR

U00 = Un_SUB(Un_SUB.NumField==1,:);
U01 = Un_SUB(Un_SUB.NumField>1,:);
U10 = Un_CA1(Un_CA1.NumField==1,:);
U11 = Un_CA1(Un_CA1.NumField>1,:);

figure; 
subplot(1,2,1); hold on
x00 = U00.RipPartRate_all; x01 = U01.RipPartRate_all;
x10 = U10.RipPartRate_all; x11 = U11.RipPartRate_all;

dat = [nanmean(x00) nanmean(x01) nanmean(x10) nanmean(x11)];
err = [nanstd(x00)/sqrt(sum(~isnan(x00))) nanstd(x01)/sqrt(sum(~isnan(x01))) nanstd(x10)/sqrt(sum(~isnan(x10))) nanstd(x11)/sqrt(sum(~isnan(x11)))];
bar(dat)
errorbar(dat,err)
ylim([0 .6])
title('Ripple Participation Rate, all SWRs')

subplot(1,2,2); hold on
x00 = U00.RipPartRate_NonSp; x01 = U01.RipPartRate_NonSp;
x10 = U10.RipPartRate_NonSp; x11 = U11.RipPartRate_NonSp;

dat = [nanmean(x00) nanmean(x01) nanmean(x10) nanmean(x11)];
err = [nanstd(x00)/sqrt(sum(~isnan(x00))) nanstd(x01)/sqrt(sum(~isnan(x01))) nanstd(x10)/sqrt(sum(~isnan(x10))) nanstd(x11)/sqrt(sum(~isnan(x11)))];
bar(dat)
errorbar(dat,err)
ylim([0 .6])
title('Ripple Participation Rate, task-related SWRs')
% % yticks([0:al:0.2])
[p,h,stats] = ranksum(x00,x01)
title(['p=' jjnum2str(p,3) ', t=' jjnum2str(stats.zval,3)])
%%
data = [x00;x01;x10;x11];
g1 = [repmat('SUB',[size([x00;x01],1) 1]); repmat('CA1',[size([x10;x11],1) 1])]
g2 = [repmat('Het',[size([x00],1) 1]); repmat('Hom',[size([x01],1) 1]);repmat('Het',[size([x10],1) 1]); repmat('Hom',[size([x11],1) 1])];

[p,tbl,stats] = anovan(data, {g1,g2}, 'model', 'interaction');
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



%% hist_RDI_session
    U0 = UnF_SUB;     U1 = UnF_CA1;
sess = unique([Un_CA1(:,2:3)]);
for s=1:size(sess,1)

temp_CA1 = U1(U1.rat==sess.rat(s) & U1.session==sess.session(s),:);
temp_SUB = U0(U0.rat==sess.rat(s) & U0.session==sess.session(s),:);
figure;
sgtitle([num2str(sess.rat(s)) '-' num2str(sess.session(s))])
varList = {'LScene','RScene','LR'};
for v=1:3
subplot(2,3,v); hold on
c0 = histogram((temp_SUB.(['RDI_' varList{v}])),'binwidth',.1,'Normalization','probability');
c0.FaceColor = CList(1,:);
xlim([-1.2 1.2]);ylim([0 .5])

subplot(2,3,v+3); hold on
c1 = histogram((temp_CA1.(['RDI_' varList{v}])),'binwidth',.1,'Normalization','probability');
c1.FaceColor = CList(2,:);
xlim([-1.2 1.2]);ylim([0 .5])

% [h,p,r] = kstest2((temp_CA1.(['RDI_' varList{v}])),(temp_SUB.(['RDI_' varList{v}])))
% 
% 
% title(['p = ' num2str(p)])

end

saveas(gca,['D:\HPC-SWR project\Processed Data_231113\Manuscript\supple\RDI_dist\' num2str(sess.rat(s)) '-' num2str(sess.session(s)) '.png'])
saveas(gca,['D:\HPC-SWR project\Processed Data_231113\Manuscript\supple\RDI_dist\' num2str(sess.rat(s)) '-' num2str(sess.session(s)) '.svg'])
close all
end