% UnitPropertySheet_Lv1


Initial_SWRFilter_common;
ROOT.Save = [ROOT.Save '\units_mat\ProfilingSheet\U1'];

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx']);
exper = {'JS','SEB','LSM'};
TargRegion = 'CA1';

for cid = 452:size(Cluster_List,1)
if ismember(Cluster_List.experimenter{cid},exper) 
        clusterID = Cluster_List.ID{cid};
        sid = find(Session_List.rat==Cluster_List.rat(cid) & Session_List.session==Cluster_List.session(cid));
[fig,nspks,onmazeAvgFR] = JMGetClusterQuals_4sheet(clusterID, ROOT,Session_List.type{sid},Cluster_List.experimenter{cid});
Cluster_List.onMazeAvgFR(cid) = onmazeAvgFR;
Cluster_List.nSpks(cid) = nspks;
end
end

writetable(Cluster_List_CA1,[ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx'],'WriteMode', 'replacefile')
