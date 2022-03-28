% UnitPropertySheet_Lv1


Initial_SWRFilter_common;
ROOT.Save = [ROOT.Save '\units_mat\ProfilingSheet\U1'];

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx']);
exper = {'JS','SEB','LSM'};
TargRegion = 'CA1';

for cid = 1:size(Cluster_List,1)
if ismember(Cluster_List.experimenter{cid},exper) 
        clusterID = Cluster_List.ID{cid};
        sid = find(Session_List.rat==Cluster_List.rat(cid) & Session_List.session==Cluster_List.session(cid));
[fig,nspks,onmazeAvgFR] = JMGetClusterQuals_4sheet(clusterID, ROOT,Session_List.type{sid},Cluster_List.experimenter{cid});
Cluster_List.onMazeAvgFR(cid) = onmazeAvgFR;
Cluster_List.nSpks(cid) = nspks;
end
end

writetable(Cluster_List,[ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx'],'WriteMode', 'replacefile')

Cluster_List_CA1 = Cluster_List;
% average peak-to-valley amplitude of waveforms ≥ 75uV 
id = Cluster_List_CA1.AvgPeaktoValley>=75;
Cluster_List_CA1 = Cluster_List_CA1(id,:);
 
% proportion of spikes within a 1ms refractory period < 1% (of total spikes)
id = Cluster_List_CA1.withinRef<1;
Cluster_List_CA1 = Cluster_List_CA1(id,:);
 
% average firing rate during the outbound journey on the stem and arms ≥ 1Hz
id = Cluster_List_CA1.onMazeAvgFR>=0.5;
Cluster_List_CA1 = Cluster_List_CA1(id,:);
 
% fast-spiking neurons (mean firing rate ≥ 10 Hz; width of the average waveform < 325us) were excluded
id = Cluster_List_CA1.onMazeAvgFR>=10 & Cluster_List_CA1.SpkWidth<325;
Cluster_List_CA1 = Cluster_List_CA1(~id,:);

writetable(Cluster_List_CA1,[ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx'],'WriteMode', 'replacefile')