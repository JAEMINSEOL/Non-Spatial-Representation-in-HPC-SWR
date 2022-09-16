% UnitPropertySheet_Lv1


Initial_SWRFilter_common;
ROOT.Save = [ROOT.Save '\units_mat\ProfilingSheet\U1'];

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx']);
exper = {'JS','SEB'};
TargRegion = 'CA1';

for cid = 1:size(Cluster_List,1)
if ismember(Cluster_List.experimenter{cid},exper) 
        clusterID = Cluster_List.ID{cid};
        sid = find(Session_List.rat==Cluster_List.rat(cid) & Session_List.session==Cluster_List.session(cid));
        try
[fig,nspks,onmazeAvgFR] = JMGetClusterQuals_4sheet(clusterID, ROOT,Session_List.type{sid},Cluster_List.experimenter{cid});
Cluster_List.onMazeAvgFR(cid) = onmazeAvgFR;
Cluster_List.nSpks(cid) = nspks;
        catch
            disp([clusterID ' is failed!'])
        end
end
end

% writetable(Cluster_List_CA1,[ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx'],'WriteMode', 'replacefile')
