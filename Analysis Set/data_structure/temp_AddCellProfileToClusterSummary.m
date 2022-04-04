
Initial_SWRFilter_common;


Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_CA1.xlsx']);
exper = {'LSM','JS','SEB'};
TargRegion = 'CA1';

for cid = 1:size(Cluster_List,1)
    if ismember(Cluster_List.experimenter{cid}, exper)
        try
        clusterID = Cluster_List.ID{cid};
 [nSPKS, FRRate,withinREFRACPortion, max_width, max_peak, max_amp, peak_ratio, LRATIO, ISODIST,  LogISIPEAKTIME,...
     valley_slope, valley_proportion,SpaInfoScore,onmazeAvgFR] = JMGetClusterQuals_lite(clusterID, ROOT,Cluster_List.experimenter{cid});
    
 Cluster_List.nSpks(cid) = nSPKS;
 Cluster_List.withinRef(cid) = withinREFRACPortion;
    Cluster_List.SpkWidth(cid) = max_width;
    Cluster_List.LogISIPeakTime(cid) = LogISIPEAKTIME;
    Cluster_List.AvgPeaktoValley(cid) = max_amp;
    Cluster_List.PeakfromBaseline(cid) = max_peak;
    Cluster_List.PeakRatio(cid) = peak_ratio;
    Cluster_List.AvgFR(cid) = FRRate;
    Cluster_List.onMazeAvgFR(cid) = onmazeAvgFR;
     Cluster_List.SI(cid) = SpaInfoScore;
        catch
        end
    end
end



writetable(Cluster_List,[ROOT.Info '\ClusterList_SWR_CA1.xlsx'],'WriteMode', 'replacefile')
Cluster_List_CA1 = Cluster_List;
% average peak-to-valley amplitude of waveforms ≥ 75uV 
id = Cluster_List_CA1.AvgPeaktoValley>=75;
Cluster_List_CA1 = Cluster_List_CA1(id,:);
 
% proportion of spikes within a 1ms refractory period < 1% (of total spikes)
id = Cluster_List_CA1.withinRef<1;
Cluster_List_CA1 = Cluster_List_CA1(id,:);
 
% average firing rate during the outbound journey on the stem and arms ≥ 1Hz
id = Cluster_List_CA1.AvgFR>=0.5;
Cluster_List_CA1 = Cluster_List_CA1(id,:);
 
% fast-spiking neurons (mean firing rate ≥ 10 Hz; width of the average waveform < 325us) were excluded
id = Cluster_List_CA1.AvgFR>=10 & Cluster_List_CA1.SpkWidth<325;
Cluster_List_CA1 = Cluster_List_CA1(~id,:);


writetable(Cluster_List_CA1,[ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx'],'WriteMode', 'replacefile')

