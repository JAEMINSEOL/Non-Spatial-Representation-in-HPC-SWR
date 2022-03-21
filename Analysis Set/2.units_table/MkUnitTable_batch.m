% MkUnitTable_batch

Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R1'];
ROOT.Unit = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';

UnitsTable_all =table;

fd = dir(ROOT.Rip);

for i=1:size(fd,1)
    id = fd(i).name;
    
    if length(id)==14
        thisSID = [id(1:6)];
    
        Recording_region_TT = Recording_region({thisSID},:);
        if strcmp(thisRegion, 'CA3')
            TargetTT = find(cellfun(cellfindn2(thisRegion),table2array(Recording_region_TT)'));
        else
            TargetTT = find(cellfun(cellfind(thisRegion),table2array(Recording_region_TT)'));
        end
        
        Spike = LoadSpikeData(ROOT, thisSID, TargetTT,cellfindn);
        load([ROOT.Behav '\' thisSID '.mat'])
        UnitsTable = MkUnitTable(ROOT,Behav,Spike,thisSID,thisRegion,TargetTT);
        
        save([ROOT.Unit '\' thisSID '.mat'], 'UnitsTable')
        UnitsTable_all = [UnitsTable_all; UnitsTable];
    end
end

writetable(UnitsTable_all,[ROOT.Save '\UnitsTable_' thisRegion '.xlsx']);

%% unit filtering

units = UnitsTable_all;
% average peak-to-valley amplitude of waveforms >= 75uV 
id = units.Peak_peakToValley_>=75;
units = units(id,:);

% proportion ofspikes within a 1ms refractory period < 1% of total spikes
id = units.WithinRefractory_perc_<1;
units = units(id,:);

% average firing rate during the outbound journey on the stem and arms >= 0.5Hz
id = units.onmazeAvgFR1DOut>=0.5;
units = units(id,:);

% fast-spiking neurons (mean firing rate >=10 Hz; width of the average waveform < 325us) were excluded
id = units.MeanFR>=10 & units.SpkWidth_maxAmp_<325;
units = units(~id,:);

UnitsTable_filtered = units;
writetable(UnitsTable_filtered,[ROOT.Save '\UnitsTable_filtered2_' thisRegion '.xlsx']);
