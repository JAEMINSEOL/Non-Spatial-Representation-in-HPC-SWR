Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';



fd = dir(ROOT.Old);

RipplesTable_all = table;

for i=1:size(fd,1)
    id = fd(i).name;
    
    if length(id)>5
        load([ROOT.Old '\' id]);
        
        thisSID.full = [id(1:6)];
        
        Recording_region_TT = Recording_region({thisSID.full},:);
        if strcmp(thisRegion, 'CA3')
            TargetTT = find(cellfun(cellfindn2(thisRegion),table2array(Recording_region_TT)'));
        else
            TargetTT = find(cellfun(cellfind(thisRegion),table2array(Recording_region_TT)'));
        end
        
        EEG = LoadEEGData(ROOT, thisSID.full, TargetTT,Params,Params_Ripple);
        
        Rip=table;
        Rip.RippleID = RipplesTable.RippleID;
        Rip.Region =  RipplesTable.Region;
        Rip.StartTime =  RipplesTable.RippleStart;
        Rip.EndTime = RipplesTable.RippleEnd;
        Rip.Duration = RipplesTable.RippleDuration;
        Rip.Gap = RipplesTable.RippleGap;
        Rip.NumAllTT = RipplesTable.NumAllTT;
        Rip.NumTT = RipplesTable.NumTT;
        
        
        
        for j=1:size(Rip,1)
            thisRippleTable = RipplesTable(j,:);
            Idx = [thisRippleTable.RippleStartIndex:thisRippleTable.RippleEndIndex];
            
            EEG_Prop=struct;
            for t=1:length(TargetTT)
                thisEEG = EEG.(['TT' num2str(TargetTT(t))]);
                [EEG_Prop.Amp(t), EEG_Prop.Freq(t), EEG_Prop.power_ripple(t), EEG_Prop.power_theta(t)] = LFP_Properties(thisEEG,Idx);
            end
            Rip.MaxVoltage(j) = max(EEG_Prop.Amp);
            Rip.MeanVoltage(j) = mean(EEG_Prop.Amp);
            Rip.MaxFreq(j) = max(EEG_Prop.Freq);
            Rip.MeanFreq(j) = mean(EEG_Prop.Freq);
            Rip.RipplePower(j) = mean(EEG_Prop.power_ripple);
            Rip.ThetaPower(j) = mean(EEG_Prop.power_theta);
        end
        
        RipplesTable_all = [RipplesTable_all; Rip];
        save([ROOT.Save '\' thisSID.full '_' thisRegion '.mat'], 'Rip')
    end
    
end

save([ROOT.Save '\RipplesTable.mat'], 'RipplesTable_all')
