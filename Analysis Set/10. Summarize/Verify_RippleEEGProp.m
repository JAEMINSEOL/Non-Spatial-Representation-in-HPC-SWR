Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R1'];
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R1'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
thisRegion0 = 'CA1';
thisRegion = 'CA1';

RippleList = readtable([ROOT.Old '\RipplesTable_' thisRegion0 '.xlsx'],'ReadRowNames',false);


Experimenter = {'LSM'};
thisSID_old='';

for sid=1:size(RippleList,1)
        if ismember(RippleList.experimenter(sid),Experimenter) 
        thisRSID = [jmnum2str(RippleList.rat(sid),3) '-' jmnum2str(RippleList.session(sid),2)];
        Recording_region_TT = Recording_region({thisRSID},:);
        
        
        temp_r = decimalToBinaryVector(RippleList.TTs(sid));
        RippleTT = flip(length((temp_r))+1 - find(temp_r))';
        
        if ~strcmp(thisRSID, thisSID_old)
            TargetTT = GetTargetTT(ROOT,thisRSID,thisRegion,Params,1);
            disp([thisRSID ' EEG data loading...'])
            EEG = LoadEEGData(ROOT, thisRSID, TargetTT,Params,Params_Ripple);
            thisSID_old = thisRSID;
        end
        
        EEG_Prop=struct;
        Idx = [RippleList.STindex(sid):RippleList.EDindex(sid)];
        Idx_pre = [RippleList.STindex(sid)-0.02*2000:RippleList.STindex(sid)];
        Idx_out = [RippleList.STindex(sid)-0.02*2000:RippleList.STindex(sid),...
            RippleList.EDindex(sid):RippleList.EDindex(sid)+0.02*2000];

        for t=1:length(TargetTT)
            try
            thisEEG = EEG.(['TT' num2str(TargetTT(t))]);
            [EEG_Prop.Amp(t), EEG_Prop.Amp_f(t), EEG_Prop.Freq(t), EEG_Prop.power_ripple(t), EEG_Prop.power_theta(t)] = LFP_Properties(thisEEG,Idx);
             [OEEG_Prop.Amp(t),OEEG_Prop.Amp_f(t) , ~, ~, ~] = LFP_Properties(thisEEG,Idx_out);
 [REEG_Prop.Amp(t),PEEG_Prop.Amp_f(t) , ~, ~, ~] = LFP_Properties(thisEEG,Idx_pre);

            catch
                 EEG_Prop.Amp(t) = nan; EEG_Prop.Amp_f(t) = nan; EEG_Prop.Freq(t)= nan; EEG_Prop.power_ripple(t)= nan; EEG_Prop.power_theta(t)= nan;
            end
        end

RippleList.MeanRaw_out(sid) = nanmean(OEEG_Prop.Amp);
RippleList.MeanFilteredPeak_out(sid) = nanmean(OEEG_Prop.Amp_f);


        end
end

RippleList.NormFilteredPeak = RippleList.MeanFilteredPeak ./ RippleList.MeanFilteredPeak_out;
RippleList.NormFilteredPeak_pre = RippleList.MeanFilteredPeak ./ nanmean(PEEG_Prop.Amp_f);
writetable(RippleList,[ROOT.Save '\RipplesTable_' thisRegion0  '.xlsx'],'WriteMode', 'replacefile')

histogram(RippleList.NormFilteredPeak)
line([1.2 1.2],[0 1600],'color','r')