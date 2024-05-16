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

%% line_RPR
U0 = Un_SUB;
U1 = Un_CA1;

U0S = U0(U0.NumField==1,:);
U0M = U0(U0.NumField>1,:);

U1S = U1(U1.NumField==1,:);
U1M = U1(U1.NumField>1,:);

x0_AS = U0S.RipPartRate_all; x0_AM = U0M.RipPartRate_all; x0_NS = U0S.RipPartRate_NonSp; x0_NM = U0M.RipPartRate_NonSp;
x1_AS = U1S.RipPartRate_all; x1_AM = U1M.RipPartRate_all; x1_NS = U1S.RipPartRate_NonSp; x1_NM = U1M.RipPartRate_NonSp;

[p,~,stats] = ranksum(x1_AS,x1_NS)

figure;
subplot(1,2,1); hold on

dat = [nanmean(x0_AS) nanmean(x0_NS)];
err = [nanstd(x0_AS)/sqrt(sum(~isnan(x0_AS))) nanstd(x0_NS)/sqrt(sum(~isnan(x0_NS))) ];
errorbar(dat,err,'color','k')
plot(dat,'linewidth',2,'color',hex2rgb('d86d35'))

dat = [nanmean(x0_AM) nanmean(x0_NM)];
err = [nanstd(x0_AM)/sqrt(sum(~isnan(x0_AM))) nanstd(x0_NM)/sqrt(sum(~isnan(x0_NM))) ];
errorbar(dat,err,'color','k')
plot(dat,'linewidth',2,'color',hex2rgb('7e4b8e'))

ylim([0 1])
title('Ripple Participation Rate, SUB')
xlim([0.5 2.5]); xticks([1:2]); xticklabels({'All','Task-related'})
ylabel('Ripple Part. Rate')

subplot(1,2,2); hold on

dat = [nanmean(x1_AS) nanmean(x1_NS)];
err = [nanstd(x1_AS)/sqrt(sum(~isnan(x1_AS))) nanstd(x1_NS)/sqrt(sum(~isnan(x1_NS))) ];
errorbar(dat,err,'color','k')
plot(dat,'linewidth',2,'color',hex2rgb('d86d35'))

dat = [nanmean(x1_AM) nanmean(x1_NM)];
err = [nanstd(x1_AM)/sqrt(sum(~isnan(x1_AM))) nanstd(x1_NM)/sqrt(sum(~isnan(x1_NM))) ];
errorbar(dat,err,'color','k')
plot(dat,'linewidth',2,'color',hex2rgb('7e4b8e'))

ylim([0 1])
title('Ripple Participation Rate, CA1')
xlim([0.5 2.5]); xticks([1:2]); xticklabels({'All','Task-related'})
ylabel('Ripple Part. Rate')
%% scatter_RPR_SI_nonSp
U0 = Un_SUB;
U1 = Un_CA1;

U0S = U0( U0.NumField==1,:);
U0M = U0( U0.NumField>1,:);

U1S = U1(U1.NumField==1,:);
U1M = U1(U1.NumField>1,:);

x0_S = U0S.RipPartRate_NonSp; x0_M = U0M.RipPartRate_NonSp;
y0_SL = U0S.RDI_LScene; y0_SR = U0S.RDI_LScene; y0_SC = U0S.RDI_LR;
y0_ML = U0M.RDI_LScene; y0_MR = U0M.RDI_LScene; y0_MC = U0M.RDI_LR;

x1_S = U1S.RipPartRate_NonSp; x1_M = U1M.RipPartRate_NonSp;
y1_SL = U1S.RDI_LScene; y1_SR = U1S.RDI_LScene; y1_SC = U1S.RDI_LR;
y1_ML = U1M.RDI_LScene; y1_MR = U1M.RDI_LScene; y1_MC = U1M.RDI_LR;

figure('position',[148,137,1208,715]);  
sgtitle('Left Scene Selectivity vs. SWR Part. Rate (task-related reactivation)')
subplot(1,2,1); hold on
scatter(y0_SL,x0_S,40,hex2rgb('d86d35'))
scatter(y0_ML,x0_M,40,hex2rgb('7e4b8e'))
xlim([-1.5 1.5]); ylim([0 1]); xlabel('selectivity index'); ylabel('SWR Part. Rate')
title('SUB'); legend({'SF','MF'})

subplot(1,2,2); hold on
scatter(y1_SL,x1_S,40,hex2rgb('d86d35'))
scatter(y1_ML,x1_M,40,hex2rgb('7e4b8e'))
xlim([-1.5 1.5]); ylim([0 1]); xlabel('selectivity index'); ylabel('SWR Part. Rate')
title('CA1'); legend({'SF','MF'})


figure('position',[148,137,1208,715]); 
sgtitle('Right Scene Selectivity vs. SWR Part. Rate (task-related reactivation)')
subplot(1,2,1); hold on
scatter(y0_SR,x0_S,40,hex2rgb('d86d35'))
scatter(y0_MR,x0_M,40,hex2rgb('7e4b8e'))
xlim([-1.5 1.5]); ylim([0 1]); xlabel('selectivity index'); ylabel('SWR Part. Rate')
title('SUB'); legend({'SF','MF'})

subplot(1,2,2); hold on
scatter(y1_SR,x1_S,40,hex2rgb('d86d35'))
scatter(y1_MR,x1_M,40,hex2rgb('7e4b8e'))
xlim([-1.5 1.5]); ylim([0 1]); xlabel('selectivity index'); ylabel('SWR Part. Rate')
title('CA1'); legend({'SF','MF'})


figure('position',[148,137,1208,715]); 
sgtitle('Choice Selectivity vs. SWR Part. Rate (task-related reactivation)')
subplot(1,2,1); hold on
scatter(y0_SC,x0_S,40,hex2rgb('d86d35'))
scatter(y0_MC,x0_M,40,hex2rgb('7e4b8e'))
xlim([-1.5 1.5]); ylim([0 1]); xlabel('selectivity index'); ylabel('SWR Part. Rate')
title('SUB'); legend({'SF','MF'})

subplot(1,2,2); hold on
scatter(y1_SC,x1_S,40,hex2rgb('d86d35'))
scatter(y1_MC,x1_M,40,hex2rgb('7e4b8e'))
xlim([-1.5 1.5]); ylim([0 1]); xlabel('selectivity index'); ylabel('SWR Part. Rate')
title('CA1'); legend({'SF','MF'})

%%
[~,ia,ic] = unique(UnitPair_field.CA1.UID1);

UF1 = UnitPair_field.CA1(ia,:);

[~,ia,ic] = unique(UnitPair_field.SUB.UID1);

UF0 = UnitPair_field.SUB(ia,:);


    UF0.Nsp_M = nanmax(abs([UF0.Nsp_L,UF0.Nsp_R,UF0.Nsp_C]),[],2);
    UF1.Nsp_M = nanmax(abs([UF1.Nsp_L,UF1.Nsp_R,UF1.Nsp_C]),[],2);
 
    Reg='CA1';
for u=1:size(UnitsTable_field.(Reg),1)
    idx = find(strcmp(UnitsTable_field.(Reg).ID(u),UnitPair_field.(Reg).UID1),1);
    if ~isempty(idx)
    UnitsTable_field.(Reg).RPR(u) = UnitPair_field.(Reg).p(idx) / UnitPair_field.(Reg).p0(idx);
    else
        UnitsTable_field.(Reg).RPR(u)=nan;
    end
end
%%
UF1 = UnitsTable_field.CA1;
UF0 = UnitsTable_field.SUB;

figure;

subplot(1,2,1)
x=abs(UF0.RDI_RScene); y=UF0.RPR;
x(~(y<=1))=[]; y(~(y<=1))=[];
y(isnan(x))=[]; x(isnan(x))=[];
scatter(x,y,40,CList(1,:),'filled')
ylim([0 1]); xlim([0 1.2])
X = [ones(size(x)) x];
B = X\y;
Rsq = 1 - sum((y - X*B).^2)/sum((y - mean(y)).^2);
title(num2str(Rsq))

subplot(1,2,2)
x=abs(UF1.RDI_RScene); y=UF1.RPR;
x(~(y<=1))=[]; y(~(y<=1))=[];
y(isnan(x))=[]; x(isnan(x))=[];
scatter(x,y,40,CList(2,:),'filled')
ylim([0 1]); xlim([0 1.2])


corrcoef(x,y)
X = [ones(size(x)) x];
B = X\y;
Rsq = 1 - sum((y - X*B).^2)/sum((y - mean(y)).^2);
title(num2str(Rsq))