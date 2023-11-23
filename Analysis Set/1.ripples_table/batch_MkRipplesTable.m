Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Processed '\ripples_mat\R0'];
ROOT.Save = [ROOT.Processed '\ripples_mat\R1'];
ROOT.R2 = [ROOT.Processed '\ripples_mat\R2'];

% ROOT.Old = [ROOT.Processed];
% ROOT.Save = [ROOT.Processed];
% ROOT.R2 = [ROOT.Processed '\ripples_mat\R2'];
% 
% ROOT.Old = [ROOT.R2];
% ROOT.Save = [ROOT.R2];


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
RegionList = {'CA1','SUB'};

for regID = 1:2
thisRegion0 = RegionList{regID};
thisRegion = RegionList{regID};

RippleList_old = readtable([ROOT.Old '\RipplesList_' thisRegion0 '_overlap.xlsx'],'ReadRowNames',false);
% RippleList_old = readtable([ROOT.Old '\RipplesTable_Behav_' thisRegion0 '.xlsx'],'ReadRowNames',false);
% RippleList_old = readtable([ROOT.Old '\RipplesTable_' thisRegion0 '_forAnalysis.xlsx'],'ReadRowNames',false);


Experimenter = {'SEB','LSM','JS'};
Experimenter = {'LSM'};
RippleList = RippleList_old;

fd = dir(ROOT.Old);

RipplesTable_all = table;
thisSID_old = '';
if ~exist(ROOT.Save), mkdir(ROOT.Save); end
RippleList.RippleDuration=(RippleList.EDtime-RippleList.STtime);
RippleList= RippleList(strcmp(RippleList.experimenter,'LSM'),:);
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
        
        for t=1:length(TargetTT)
            try
            thisEEG = EEG.(['TT' num2str(TargetTT(t))]);
            [EEG_Prop.Amp(t), EEG_Prop.Amp_f(t), EEG_Prop.Freq(t), EEG_Prop.power_ripple(t), EEG_Prop.power_theta(t)] = LFP_Properties(thisEEG,Idx);
            catch
                 EEG_Prop.Amp(t) = nan; EEG_Prop.Amp_f(t) = nan; EEG_Prop.Freq(t)= nan; EEG_Prop.power_ripple(t)= nan; EEG_Prop.power_theta(t)= nan;
            end
        end
        RippleList.MeanFilteredPeak_all(sid) = nanmean(EEG_Prop.Amp_f);
        RippleList.StdFilteredPeak_all(sid) = nanstd(EEG_Prop.Amp_f);
        
%         EEG_Prop=struct;
%         Idx = [RippleList.STindex(sid):RippleList.EDindex(sid)];
%         
%         for t=1:length(RippleTT)
%             try
%             thisEEG = EEG.(['TT' num2str(RippleTT(t))]);
%             [EEG_Prop.Amp(t), EEG_Prop.Amp_f(t), EEG_Prop.Freq(t), EEG_Prop.power_ripple(t), EEG_Prop.power_theta(t)] = LFP_Properties(thisEEG,Idx);
%             catch
%                  EEG_Prop.Amp(t) = nan; EEG_Prop.Amp_f(t) = nan; EEG_Prop.Freq(t)= nan; EEG_Prop.power_ripple(t)= nan; EEG_Prop.power_theta(t)= nan;
%             end
%         end
        RippleList.MeanFilteredPeak(sid) = nanmean(EEG_Prop.Amp_f);
        RippleList.StdFilteredPeak(sid) = nanstd(EEG_Prop.Amp_f);
        RippleList.MeanRaw(sid) = nanmean(EEG_Prop.Amp);
        RippleList.MaxFreq(sid) = nanmax(EEG_Prop.Freq);
        RippleList.MeanFreq(sid) = nanmean(EEG_Prop.Freq);
        RippleList.RipplePower(sid) = nanmean(EEG_Prop.power_ripple);
        RippleList.ThetaPower(sid) = nanmean(EEG_Prop.power_theta);
%         writetable(Rip,[ROOT.Save '\RipplesTable_' thisSID '_' thisRegion '.xlsx'], 'overwritesheet')
    end
    
end
RippleList(RippleList.MeanFilteredPeak_all==0,:)=[];
RippleList(RippleList.rat==427,:)=[];
writetable(RippleList,[ROOT.Save '\RipplesTable_' thisRegion  '_forAnalysis.xlsx'],'WriteMode', 'replacefile')
%  load([ROOT.Save '\RipplesTable.mat'])
% writetable(RippleList,[ROOT.Save  '\RipplesTable_Behav_' thisRegion0 '.xlsx'],'WriteMode', 'replacefile')
end

