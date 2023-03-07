% MkUnitTable_batch

Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data (Blur)'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R1'];
ROOT.Unit = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';

UnitsTable_all =table;

Cluster_List = readtable([ROOT.Save '\units_mat\U0\ClusterList_SWR_CA1_filtered.xlsx']);
exper = {'JS','LSM','SEB'};
exper2 = {'JS','LSM'};
% exper = {'SEB'};
TargRegion = 'CA1';


for cid=1:size(Cluster_List,1)
    id = Cluster_List.ID{cid};
    if ~ismember(Cluster_List.experimenter{cid},exper), continue; end
 
        thisSID = id(1:6);
    
        Recording_region_TT = Recording_region({thisSID},:);
        if strcmp(thisRegion, 'CA3')
            TargetTT = find(cellfun(Params.cellfindn2(thisRegion),table2array(Recording_region_TT)'));
        else
            TargetTT = find(cellfun(Params.cellfindn2(thisRegion),table2array(Recording_region_TT)'));
        end
        
        Spike = LoadSpikeData(ROOT, thisSID, TargetTT,Params.cellfindn);
        load([ROOT.Behav '\' thisSID '.mat'])
        
        UnitsTable = MkUnitTable_blur(ROOT,Behav,Spike,id,thisRegion,TargetTT,Params);
        
        save([ROOT.Unit '\' thisSID '.mat'], 'UnitsTable')
        UnitsTable_all = [UnitsTable_all; UnitsTable];
        disp([id ' is finished!'])
end

writetable(UnitsTable_all,[ROOT.Unit '\UnitsTable3_' thisRegion '.xlsx'],'writemode','overwrite');

%% unit filtering

units = UnitsTable_all;
% average peak-to-valley amplitude of waveforms >= 75uV 
id = units.AvgPeaktoValley>=75;
units = units(id,:);

% proportion ofspikes within a 1ms refractory period < 1% of total spikes
id = units.withinRef<1;
units = units(id,:);

% fast-spiking neurons (mean firing rate >=10 Hz; width of the average waveform < 325us) were excluded
id = units.AvgFR>=10 & units.SpkWidth<325;
units = units(~id,:);

UnitsTable_filtered = units;
writetable(UnitsTable_filtered,[ROOT.Unit '\UnitsTable_filtered_' thisRegion '.xlsx'],'writemode','overwrite');

% average firing rate during the outbound journey on the stem and arms >= 0.5Hz
id = units.AvgFR>=0.5;
units = units(id,:);

id = ismember(units.experimenter,exper2);
units = units(id,:);

id = nanmax(abs([units.ReMap_No units.ReMap_Lo units.ReMap_Hi]),[],2)>.5;
units = units(id,:);

UnitsTable_anal = units;
writetable(UnitsTable_anal,[ROOT.Unit '\UnitsTable_' thisRegion '_forAnalysis.xlsx'],'writemode','overwrite');
