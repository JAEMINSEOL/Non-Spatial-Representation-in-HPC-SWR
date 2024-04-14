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

%%
UnitsTable_temp.SUB_SF = Un_SUB(Un_SUB.NumField==1,:);
UnitsTable_temp.SUB_MF = Un_SUB(Un_SUB.NumField>1,:);
UnitsTable_temp.CA1_SF = Un_CA1(Un_CA1.NumField==1,:);
UnitsTable_temp.CA1_MF = Un_CA1(Un_CA1.NumField>1,:);
%% for fig 2
crit_list = [0.05:0.05:0.4];
unit_list = {'SUB_SF','SUB_MF','CA1_SF','CA1_MF'};

figure('Position',[86,165,1215,813])
for ulid = 1:4
    subplot(2,2,ulid); hold on
    Un_temp = UnitsTable_temp.(unit_list{ulid});
for cri = 1:length(crit_list)
    crit =crit_list(cri);
    sum_all = size(Un_temp,1);
    sum_ssi(cri) = sum((abs(Un_temp.RDI_LScene)>crit | abs(Un_temp.RDI_RScene)>crit) & ~(abs(Un_temp.RDI_LR)>crit))/sum_all;
    sum_csi(cri) = sum(~(abs(Un_temp.RDI_LScene)>crit | abs(Un_temp.RDI_RScene)>crit) & (abs(Un_temp.RDI_LR)>crit))/sum_all;
    sum_cssi(cri) = sum((abs(Un_temp.RDI_LScene)>crit | abs(Un_temp.RDI_RScene)>crit) & (abs(Un_temp.RDI_LR)>crit))/sum_all;
end

plot(sum_ssi,'color',hex2rgb('e8aa42'),'Marker','o')
plot(sum_csi,'color',hex2rgb('d83f31'),'Marker','o')
plot(sum_cssi,'color',hex2rgb('025464'),'Marker','o')

xticks([1:length(crit_list)])
xticklabels(crit_list)
xlabel('Selective index criterion')
ylim([0 1])
ylabel('Selective cell proportion')

legend({'Scene only','Choice only','Scene and Choice'},'location','northoutside','Orientation','horizontal')

title(unit_list{ulid})
end

%% additional analysis
crit_list = [0.05:0.05:0.4];
unit_list = {'SUB_SF','SUB_MF','CA1_SF','CA1_MF'};

figure('Position',[86,165,1215,813])
for ulid = 1:4
    subplot(2,2,ulid); hold on
    Un_temp = UnitsTable_temp.(unit_list{ulid});
for cri = 1:length(crit_list)
    crit =crit_list(cri);
    sum_all = size(Un_temp,1);
    sum_ssi(cri) = sum((abs(Un_temp.RDI_LScene)>crit | abs(Un_temp.RDI_RScene)>crit))/sum_all;
    sum_csi(cri) = sum((abs(Un_temp.RDI_LR)>crit))/sum_all;
end

plot(sum_ssi,'color',hex2rgb('e8aa42'),'Marker','o')
plot(sum_csi,'color',hex2rgb('d83f31'),'Marker','o')

xticks([1:length(crit_list)])
xticklabels(crit_list)
xlabel('Selective index criterion')
ylim([0 1])
ylabel('Selective cell proportion')

legend({'Scene-selective','Choice-selective'},'location','northoutside','Orientation','horizontal')

title(unit_list{ulid})
end

    
