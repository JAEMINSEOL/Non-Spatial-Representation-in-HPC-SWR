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


FRMaps_SUB = LoadFRMap(ROOT,Un_SUB);
FRMaps_CA1 = LoadFRMap(ROOT,Un_CA1);
%%
crit_pf = .33;

for sid=1:size(Un_SUB,1)
[field_count, start_index, end_index, field_size, h] = getFRfields_v2_jm(FRMaps_SUB(1,:,sid),crit_pf);
Un_SUB.NumField_FR(sid) = field_count;
end

for sid=1:size(Un_CA1,1)
[field_count, start_index, end_index, field_size, h] = getFRfields_v2_jm(FRMaps_CA1(1,:,sid),crit_pf);
Un_CA1.NumField_FR(sid) = field_count;
end
%%
%% bar_RPR_het
U0 = Un_SUB;
U1 = Un_CA1;

U0_SS = U0(U0.NumField_FR==1 & U0.NumField==1,:);
U0_SM = U0(U0.NumField_FR==1 & U0.NumField>1,:);
U0_MM = U0(U0.NumField_FR>1 & U0.NumField>1,:);
U0_MS = U0(U0.NumField_FR>1 & U0.NumField<=1,:);

U1_SS = U1(U1.NumField_FR==1 & U1.NumField==1,:);
U1_SM = U1(U1.NumField_FR==1 & U1.NumField>1,:);
U1_MM = U1(U1.NumField_FR>1 & U1.NumField>1,:);
U1_MS = U1(U1.NumField_FR>1 & U1.NumField<=1,:);

figure; 
subplot(2,2,1); hold on
x0_SS = U0_SS.RipPartRate_all; x0_SM = U0_SM.RipPartRate_all; x0_MM = U0_MM.RipPartRate_all;
x1_SS = U1_SS.RipPartRate_all; x1_SM = U1_SM.RipPartRate_all; x1_MM = U1_MM.RipPartRate_all;

dat = [nanmean(x0_SS) nanmean(x0_SM) nanmean(x0_MM)];
err = [nanstd(x0_SS)/sqrt(sum(~isnan(x0_SS))) nanstd(x0_SM)/sqrt(sum(~isnan(x0_SM))) nanstd(x0_MM)/sqrt(sum(~isnan(x0_MM)))];

bar(dat)
errorbar(dat,err)
ylim([0 .6])
title('Ripple Participation Rate, all SWRs - SUB')
xticks([1:3]); xticklabels({'Single-Single','Single-Multi','Multi-Multi'})
ylabel('Ripple Part. Rate')

subplot(2,2,2); hold on
dat = [nanmean(x1_SS) nanmean(x1_SM) nanmean(x1_MM)];
err = [nanstd(x1_SS)/sqrt(sum(~isnan(x1_SS))) nanstd(x1_SM)/sqrt(sum(~isnan(x1_SM))) nanstd(x1_MM)/sqrt(sum(~isnan(x1_MM)))];

bar(dat)
errorbar(dat,err)
ylim([0 .6])
title('Ripple Participation Rate, all SWRs - CA1')
xticks([1:3]); xticklabels({'Single-Single','Single-Multi','Multi-Multi'})
ylabel('Ripple Part. Rate')


subplot(2,2,3); hold on
x0_SS = U0_SS.RipPartRate_NonSp; x0_SM = U0_SM.RipPartRate_NonSp; x0_MM = U0_MM.RipPartRate_NonSp;
x1_SS = U1_SS.RipPartRate_NonSp; x1_SM = U1_SM.RipPartRate_NonSp; x1_MM = U1_MM.RipPartRate_NonSp;

dat = [nanmean(x0_SS) nanmean(x0_SM) nanmean(x0_MM)];
err = [nanstd(x0_SS)/sqrt(sum(~isnan(x0_SS))) nanstd(x0_SM)/sqrt(sum(~isnan(x0_SM))) nanstd(x0_MM)/sqrt(sum(~isnan(x0_MM)))];

bar(dat)
errorbar(dat,err)
ylim([0 .6])
title('Ripple Participation Rate, task-related SWRs - SUB')
xticks([1:3]); xticklabels({'Single-Single','Single-Multi','Multi-Multi'})
ylabel('Ripple Part. Rate')

subplot(2,2,4); hold on
dat = [nanmean(x1_SS) nanmean(x1_SM) nanmean(x1_MM)];
err = [nanstd(x1_SS)/sqrt(sum(~isnan(x1_SS))) nanstd(x1_SM)/sqrt(sum(~isnan(x1_SM))) nanstd(x1_MM)/sqrt(sum(~isnan(x1_MM)))];

bar(dat)
errorbar(dat,err)
ylim([0 .6])
title('Ripple Participation Rate, task-related SWRs - CA1')
xticks([1:3]); xticklabels({'Single-Single','Single-Multi','Multi-Multi'})
ylabel('Ripple Part. Rate')