Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R5'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R4'];
ROOT.Units = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];
if ~exist(ROOT.Fig3), mkdir(ROOT.Fig3); end
if ~exist(ROOT.Rip), mkdir(ROOT.Rip); end

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'SUB';
thisRegion2 = 'SUB';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion2 '.xlsx']);
% UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion2 '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']); 
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
randN=5000;
%%
% RipplesTable_p = RipplesTable(RipplesTable.nPCs>2,:);
RipplesTable_p = RipplesTable;
for r=1:size(RipplesTable_p,1)
    RipID = RipplesTable_p.ID{r};
    thisReact = ReactTable(strcmp(ReactTable.RippleID,RipID),:);
    thisUnits=table;
    for u=1:size(thisReact,1)
        thisUnits =[thisUnits;UnitsTable_A(find(cellfun(Params.cellfind(thisReact.UnitID(u)),UnitsTable_A.ID)),:)];
    end

    thisUnits_L = thisUnits(~isnan(thisUnits.RDI_LScene),:);
    thisUnits_R = thisUnits(~isnan(thisUnits.RDI_RScene),:);
    thisUnits_C = thisUnits(~isnan(thisUnits.RDI_LR),:);

    

    RipplesTable_p.nRDI_MF(r) = nanmean(thisUnits.NumField>1);
    RipplesTable_p.nRDI_hetero_SC(r) = nanmean(thisUnits.MultiVar_SC);
    RipplesTable_p.mRDI_hetero_SC(r) = nanmean(thisUnits.RDI_hetero_SC);

    RipplesTable_p.nRDI_hetero_L(r) = nanmean(thisUnits.MultiVar_L);
    RipplesTable_p.mRDI_hetero_L(r) = nanmean(thisUnits.RDI_hetero_L);

    RipplesTable_p.nRDI_hetero_R(r) = nanmean(thisUnits.MultiVar_R);
    RipplesTable_p.mRDI_hetero_R(r) = nanmean(thisUnits.RDI_hetero_R);
    
    RipplesTable_p.nRDI_hetero_C(r) = nanmean(thisUnits.MultiVar_C);
    RipplesTable_p.mRDI_hetero_C(r) = nanmean(thisUnits.RDI_hetero_C);


   disp([RipID ' is finished!']) 
end
 suff = '_RDIs';
writetable(RipplesTable_p,[ROOT.Save '\RipplesTable_' thisRegion2 suff '.xlsx'],'writemode','replacefile')
