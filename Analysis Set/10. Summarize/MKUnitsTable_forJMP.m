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
%% ANOVA_RPR_SI_SFMF_FRTP_table
U0 = Un_SUB;
U1 = Un_CA1;

U0F = Un_SUB_FR;
U1F = Un_CA1_FR;

US0 = U0(U0.NumField==1,:);
US1 = U0(U0.NumField>1,:);
UC0 = U1(U1.NumField==1,:);
UC1 = U1(U1.NumField>1,:);

US0F = U0F(U0.NumField==1,:);
US1F = U0F(U0.NumField>1,:);
UC0F = U1F(U1.NumField==1,:);
UC1F = U1F(U1.NumField>1,:);

data=[];
vars='RipPartRate_all';
data.RPR_A = [US0.(vars);UC0.(vars);US1.(vars);UC1.(vars)];

vars='RipPartRate_NonSp';
data.RPR_Nsp = [US0.(vars);UC0.(vars);US1.(vars);UC1.(vars)];

vars='RDI_LScene';
data.SI_L_TP = [US0.(vars);UC0.(vars);US1.(vars);UC1.(vars)];
data.SI_L_FR = [US0F.(vars);UC0F.(vars);US1F.(vars);UC1F.(vars)];

vars='RDI_RScene';
data.SI_R_TP = [US0.(vars);UC0.(vars);US1.(vars);UC1.(vars)];
data.SI_R_FR = [US0F.(vars);UC0F.(vars);US1F.(vars);UC1F.(vars)];

vars='RDI_LR';
data.SI_C_TP = [US0.(vars);UC0.(vars);US1.(vars);UC1.(vars)];
data.SI_C_FR = [US0F.(vars);UC0F.(vars);US1F.(vars);UC1F.(vars)];

g3_s={};  shet_s=[];
for u=1:size(US1,1)
    if US1.FieldHet(u)== 0 || US1.FieldHet(u)== 1
        g3_s{u,1} = 'HOM';
        if US1.FieldHet(u)== 0
            shet_s(u,1)=0;
        else
            shet_s(u,1)=1;
        end
    elseif  US1.FieldHet(u)== 2 || US1.FieldHet(u)== 3
        g3_s{u,1} = 'HET';
        if US1.FieldHet(u)== 2
            shet_s(u,1)=0;
        else
            shet_s(u,1)=1;
        end
            else
         g3_s{u,1} = 'nan';
    end
end

g3_c={}; shet_c=[];
for u=1:size(UC1,1)
    if UC1.FieldHet(u)== 0 || UC1.FieldHet(u)== 1
        g3_c{u,1} = 'HOM';
        if UC1.FieldHet(u)== 0
            shet_c(u,1)=0;
        else
            shet_c(u,1)=1;
        end
    elseif  UC1.FieldHet(u)== 2 || UC1.FieldHet(u)== 3
        g3_c{u,1} = 'HET';
        if UC1.FieldHet(u)== 2
            shet_c(u,1)=0;
        else
            shet_c(u,1)=1;
        end
    else
         g3_c{u,1} = 'nan';
    end
end

g1 = [repmat('SUB',[size([US0],1) 1]); repmat('CA1',[size([UC0],1) 1]);...
    repmat('SUB',[size([US1],1) 1]); repmat('CA1',[size([UC1],1) 1])];
g2 = [repmat('SF',[size([US0],1) 1]); repmat('SF',[size([UC0],1) 1]); repmat('MF',[size([US1],1) 1]); repmat('MF',[size([UC1],1) 1])];

g3 = [num2cell(nan(size(US0,1),1)); num2cell(nan(size(UC0,1),1)); g3_s; g3_c];

dat_xls=table;
dat_xls.Region=g1;
dat_xls.Field=g2;
dat_xls.MFHet=g3;

dat_xls.RPR_All=data.RPR_A;
dat_xls.RPR_Nsp=data.RPR_Nsp;

dat_xls.SI_L_FR=data.SI_L_FR;
dat_xls.SI_R_FR=data.SI_R_FR;
dat_xls.SI_C_FR=data.SI_C_FR;

dat_xls.SI_L_TP=data.SI_L_TP;
dat_xls.SI_R_TP=data.SI_R_TP;
dat_xls.SI_C_TP=data.SI_C_TP;

writetable(dat_xls,[ROOT.Processed '\Units_All.xlsx'],'writemode','replacefile')