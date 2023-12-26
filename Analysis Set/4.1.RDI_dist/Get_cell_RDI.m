addpath('D:\HPC-SWR project\Analysis Program')
Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R4'];
ROOT.Units = [ROOT.Save '\units_mat\U2'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];
if ~exist(ROOT.Fig3), mkdir(ROOT.Fig3); end
if ~exist(ROOT.Rip), mkdir(ROOT.Rip); end

% parpool('local',8); % Change 4 to the number of workers you want
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


RegionList = {'CA1','SUB'};
for reg=1:2
thisR = RegionList{reg};
% 

thisRegion0 = thisR;
thisRegion = thisR;
thisRegion2 = [thisR '_field'];


% RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis' '.xlsx']);
ReactTableA = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);
ReactTableB = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion2 '.xlsx']);


UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);

TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
randN=2000;

UnitsTable = UnitsTable_B;

%%
% RipplesTable_p = RipplesTable(RipplesTable.nPCs>2,:);
RipplesTable_p = RipplesTable;
for uid=1:size(UnitsTable_A,1)
    thisUnit = UnitsTable_A(uid,:);
    thisFields = UnitsTable_B(strncmp(UnitsTable_A.ID(uid),UnitsTable_B.ID,12),:);

    if ~isempty(thisFields)
    [~,t] = max(abs(thisFields.RDI_LScene));
    UnitsTable_A.RDI_LScene(uid) = thisFields.RDI_LScene(t);

        [~,t] = max(abs(thisFields.RDI_RScene));
    UnitsTable_A.RDI_RScene(uid) = thisFields.RDI_RScene(t);

        [~,t] = max(abs(thisFields.RDI_LR));
    UnitsTable_A.RDI_LR(uid) = thisFields.RDI_LR(t);
    else
        UnitsTable_A.RDI_LScene(uid) = nan;
        UnitsTable_A.RDI_RScene(uid) = nan;
        UnitsTable_A.RDI_LR(uid) = nan;
    end
end
writetable(UnitsTable_A,[ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx'],'WriteMode','replacefile')
end
