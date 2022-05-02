Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R1'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
RippleList_old = readtable([ROOT.Old '\RipplesList_CA1_filtered.xlsx'],'ReadRowNames',false);
thisRegion = 'CA1';
Experimenter = {'LSM','JS','SEB'};

RippleList = RippleList_old;

fd = dir(ROOT.Old);

RipplesTable_all = table;
thisSID_old = '';

RippleList.RippleDuration=(RippleList.EDtime-RippleList.STtime);
for sid=1:size(RippleList,1)
    
    if ismember(RippleList.experimenter(sid),Experimenter)
        thisSID = [jmnum2str(RippleList.rat(sid),3) '-' jmnum2str(RippleList.session(sid),2)];
        Recording_region_TT = Recording_region({thisSID},:);
        TargetTT = find(cellfun(cellfind(thisRegion),table2array(Recording_region_TT)'));
        
        if ~strcmp(thisSID, thisSID_old)
            disp([thisSID ' EEG data loading...'])
            EEG = LoadEEGData(ROOT, thisSID, TargetTT,Params,Params_Ripple);
            thisSID_old = thisSID;
        end
        
        EEG_Prop=struct;
        Idx = [RippleList.STindex(sid):RippleList.EDindex(sid)];
        for t=1:length(TargetTT)
            try
            thisEEG = EEG.(['TT' num2str(TargetTT(t))]);
            [EEG_Prop.Amp(t), EEG_Prop.Freq(t), EEG_Prop.power_ripple(t), EEG_Prop.power_theta(t)] = LFP_Properties(thisEEG,Idx);
            catch
                 EEG_Prop.Amp(t) = nan; EEG_Prop.Freq(t)= nan; EEG_Prop.power_ripple(t)= nan; EEG_Prop.power_theta(t)= nan;
            end
        end
        RippleList.MaxVoltage(sid) = nanmax(EEG_Prop.Amp);
        RippleList.MeanVoltage(sid) = nanmean(EEG_Prop.Amp);
        RippleList.MaxFreq(sid) = nanmax(EEG_Prop.Freq);
        RippleList.MeanFreq(sid) = nanmean(EEG_Prop.Freq);
        RippleList.RipplePower(sid) = nanmean(EEG_Prop.power_ripple);
        RippleList.ThetaPower(sid) = nanmean(EEG_Prop.power_theta);
%         writetable(Rip,[ROOT.Save '\RipplesTable_' thisSID '_' thisRegion '.xlsx'], 'overwritesheet')
    end
    
end
writetable(RippleList,[ROOT.Save '\RipplesTable_' thisRegion '.xlsx'],'WriteMode', 'overwrite')
%  load([ROOT.Save '\RipplesTable.mat'])


