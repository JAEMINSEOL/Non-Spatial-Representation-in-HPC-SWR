Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);

Cluster_List = readtable([ROOT.Info '\ClusterSummary.xlsx']);
exper = {'SEB'};
TargRegion = 'CA1';

for cid = 1:size(Cluster_List,1)
    if ismember(Cluster_List.experimenter{cid},exper) && Cluster_List.rat(cid)>=52
        try
        clusterID = Cluster_List.ID{cid};
 [nSPKS, withinREFRACPortion, max_width, max_peak, max_amp, peak_ratio, LRATIO, ISODIST,  LogISIPEAKTIME, valley_slope, valley_proportion] = JMGetClusterQuals_lite(clusterID, ROOT);
    Cluster_List.withinRef(cid) = withinREFRACPortion;
    Cluster_List.SpkWidth(cid) = max_width;
    Cluster_List.LogISIPeakTime(cid) = LogISIPEAKTIME;
    Cluster_List.avgPeaktoValley(cid) = max_amp;
    Cluster_List.PeakfromBaseline(cid) = max_peak;
    Cluster_List.PeakRatio(cid) = peak_ratio;
        catch
        end
    end
end

writetable(Cluster_List,[ROOT.Info '\ClusterSummary.xlsx'])


Cluster_List_CA1 = readtable([ROOT.Info '\ClusterList_SWR_CA1.xlsx']);
for cid = 1:size(Cluster_List_CA1,1)
    thisCLID = Cluster_List_CA1.ID{cid};
    cid2 = find(strcmp(thisCLID,Cluster_List.ID));
    Cluster_List_CA1.nSpks(cid) = Cluster_List.nSpks(cid2);
    Cluster_List_CA1.AvgPeaktoValley(cid) = Cluster_List.avgPeaktoValley(cid2);
    Cluster_List_CA1.PeakfromBaseline(cid) = Cluster_List.PeakfromBaseline(cid2);
    Cluster_List_CA1.PeakRatio(cid) = Cluster_List.PeakRatio(cid2);
    Cluster_List_CA1.SpkWidth(cid) = Cluster_List.SpkWidth(cid2);
    Cluster_List_CA1.withinRef(cid) = Cluster_List.withinRef(cid2);
    Cluster_List_CA1.LogISIPeakTime(cid) = Cluster_List.LogISIPeakTime(cid2);
    Cluster_List_CA1.onMazeAvgFR(cid) = Cluster_List.onMazeAvgFR(cid2);
end
writetable(Cluster_List_CA1,[ROOT.Info '\ClusterList_SWR_CA1.xlsx'])

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


writetable(Cluster_List_CA1,[ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx'])
