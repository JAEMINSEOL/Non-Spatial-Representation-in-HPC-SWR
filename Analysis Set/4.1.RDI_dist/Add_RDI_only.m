% Add_RDI_only

Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed];
ROOT.Rip = [ROOT.Save '\ripples_mat\R1'];
ROOT.Unit = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1_field';

UnitsTable_all =table;

Cluster_List = readtable([ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
exper = {'LSM'};
% exper2 = {'JS','LSM'};
% exper = {'SEB'};
TargRegion = 'SUB_field';
% gcp = parpool(4);

for cid=1:size(Cluster_List,1)
    id = Cluster_List.ID{cid};
    if ~ismember(Cluster_List.experimenter{cid},exper), continue; end

    thisSID = id(1:6);

    Recording_region_TT = Recording_region({thisSID},:);

    TargetTT = find(cellfun(Params.cellfindn2(thisRegion),table2array(Recording_region_TT)'));

    if strcmp(thisRegion,'SUB_field') | strcmp(thisRegion,'CA1_field')
        Spike = LoadSpikeData_field(ROOT, thisSID, TargetTT,Params.cellfindn);
    else
        Spike = LoadSpikeData(ROOT, thisSID, TargetTT,Params.cellfindn);
    end

    
%     load([ROOT.Behav '\' thisSID '.mat'])
% 
%     UnitsTable = MkUnitTable(ROOT,Behav,Spike,id,thisRegion,TargetTT,Params);
%     save([ROOT.Unit '\' id '.mat'], 'UnitsTable')
%     UnitsTable_all = [UnitsTable_all; UnitsTable];

%%
% try
%     UnitsTable = load([ROOT.Unit '\' id '.mat']);
%     [field_count, start_index, end_index, field_size, h] = ...
%         field_boundary_function_2f2(UnitsTable.thisFieldMap1D, id);
%     UnitsTable.start_index = start_index;
%     UnitsTable.end_index = end_index;
%     UnitsTable.field_size = field_size;
%     save([ROOT.Unit '\' id '.mat'], 'UnitsTable');
% 
%     Cluster_List.FieldSize_TP(cid) = field_size;
% 
% catch
%      Cluster_List.FieldSize_TP(cid) = nan;
% end



end
%%
% writetable(Cluster_List, [ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis.xlsx'],'writemode','replacefile');
%%
% writetable(UnitsTable_all,[ROOT.Save '\UnitsTable_' thisRegion '_RDI_FR.xlsx']);

% unit = readtable([ROOT.Save '\UnitsTable_' thisRegion '_RDIperm.xlsx']);
% 
% unit = UnitsTable_all;
% 
% RDIP = min([unit.RateP_LScene, unit.RateP_RScene, unit.RateP_LR],[],2);
% sum(RDIP<0.05)
% 
% figure
% 
% hold on
% cdfplot(abs(unit.RateP_LScene))
% cdfplot(abs(unit.RateP_RScene))
% cdfplot(abs(unit.RateP_LR))
% 
% sum(abs(unit.RDI_LScene)>0.1)