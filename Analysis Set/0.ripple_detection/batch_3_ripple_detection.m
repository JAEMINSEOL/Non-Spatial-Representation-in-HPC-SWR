Initial_SWRFilter_common;

SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

Experimenter = {'LSM'};
thisRegion = 'CA1';
ROOT.Save = [ROOT.Save '\ripples_mat\R0\' thisRegion ];
if ~exist(ROOT.Save), mkdir(ROOT.Save); end

fd = dir(ROOT.Save);

RipplesTable_all = table;

for sid=1:size(SessionList,1)
    if SessionList.include(sid) && ismember(SessionList.experimenter(sid),Experimenter)
        thisRID = jmnum2str(SessionList.rat(sid),3);
        thisSID = jmnum2str(SessionList.session(sid),2);
        thisRSID = [thisRID '-' thisSID];
        
%         TargetTT = [1:24]';
        TargetTT = GetTargetTT(ROOT,thisRSID,thisRegion,Params,1);
        if isempty(TargetTT), continue;  end
        %% EEG
        EEG = LoadEEGData(ROOT, thisRSID, TargetTT,Params,Params_Ripple);
        Ripples = struct;
        for thisTTID=1:24
            try
                raw = EEG.(['TT' num2str(thisTTID)]).Raw;
                envelope_smoothed = EEG.(['TT' num2str(thisTTID)]).Envelope_smoothed;
                
                cscID = [thisRID '-' thisSID '-' num2str(thisTTID)];
                % load CSC data
                cscData = loadCSC(cscID, ROOT.Raw.Mother, Params.CSCfileTag, Params.exportMODE, Params.behExtraction,'CSC');
                Timestamps_expand = cscData.Timestamps ./ 1000000; % make to sec unit
                for i = 1 : length(Timestamps_expand)
                    for count = 2 : 512
                        Timestamps_expand(count, i) = Timestamps_expand(1, i) + (1/Params_Ripple.Fs)*(count-1);
                    end
                end
                Timestamps_expand = Timestamps_expand(:);
                %% remove noise period using session threshold
                temp_stat(1) = mean(envelope_smoothed);
                temp_stat(2) = std(envelope_smoothed,1);
                
                Params_Ripple.threshold(1) = temp_stat(1)+ Params_Ripple.noiseSTD * temp_stat(2);
                
                aboveNoiseThreshold = find(envelope_smoothed > Params_Ripple.threshold(1));
                
                envelope_smoothed_noiseRemoved = envelope_smoothed;
                envelope_smoothed_noiseRemoved(aboveNoiseThreshold,1) = NaN;
                
                %
                Params_Ripple.envelope_stat(1) = nanmean(envelope_smoothed_noiseRemoved);
                Params_Ripple.envelope_stat(2) = nanstd(envelope_smoothed_noiseRemoved,1);
                
                Params_Ripple.threshold(2) = Params_Ripple.envelope_stat(1) + Params_Ripple.thresholdSTD * Params_Ripple.envelope_stat(2);
                Params_Ripple.threshold(3) = Params_Ripple.envelope_stat(1) + Params_Ripple.beginthresholdSTD * Params_Ripple.envelope_stat(2);
                %% ripple detection
                ripples = []; ripples_index = []; % ripples_index is required for 3TToverlap
                
                beginThreshold =(envelope_smoothed > Params_Ripple.threshold(3));
                
                try
                    ripples_index = jmBin2Index(beginThreshold,1);
                    
                    
                    
                    % grouping
                    ripples2=[]; ripples_index2=[];
                    i=1;
                    while i<size(ripples_index,1)
                        j=i+1;
                        while ripples_index(j,1) - ripples_index(j-1,2) <= Params_Ripple.groupingInterval*Params_Ripple.Fs
                            j=j+1;
                            if j>size(ripples_index,1), break; end
                        end
                        ripples_index2 = [ripples_index2;[ripples_index(i,1),ripples_index(j-1,2)]];
                        i=j;
                    end
                    if ripples_index(end,2)~=ripples_index2(end,2)
                        ripples_index2 = [ripples_index2 ; ripples_index(end,1:2)];
                    end
                    ripples_index= ripples_index2;
                    
                    % ripple peak vs. ripple threshold
                    for i=1:size(ripples_index,1)
                        pk = max(envelope_smoothed(ripples_index(i,1):ripples_index(i,2))) ;
                        if pk <= Params_Ripple.threshold(2) || pk > Params_Ripple.threshold(1)
                            ripples_index(i,3)=0;
                        else
                            ripples_index(i,3)=1;
                        end
                    end
                    ripples_index(ripples_index(:,3)==0,:)=[];
                    ripples_index(:,3)=[];
                    
                    % duration check
                    ripples_index(:,3) = (ripples_index(:,2)-ripples_index(:,1));
                    ripples_index(ripples_index(:,3)<Params_Ripple.minDuration*Params.Fs,:) =[];
                    ripples_index(:,3)=[];
                    
                    % noise check
                EEGr = max(abs(raw));
                for r=1:size(ripples_index,1)
                ripples_index(r,3) = sum(abs(raw(ripples_index(r,1):ripples_index(r,2)))>=EEGr);
                end
                ripples_index(ripples_index(:,3)>1,:) =[];
                    ripples_index(:,3)=[];
                    
                    
                    ripples = zeros(2,3);
                    ripples(:,1) = Timestamps_expand(1,end);
                    ripples(:,2) = EEG.(['TT' num2str(thisTTID)]).Raw(1,end);
                    ripples(:,3) = envelope_smoothed_noiseRemoved(1,end);
                    
                    Ripples.(['TT' num2str(thisTTID)]).ripples = ripples;
                    Ripples.(['TT' num2str(thisTTID)]).index = ripples_index;
                    
                    disp([thisRID '-' thisSID '-' num2str(thisTTID) ', is finished!']);
                catch
                    disp([thisRID '-' thisSID '-' num2str(thisTTID) ', aboveThreshold is empty']);
                end
            catch
            end
        end
        save([ROOT.Save '\' thisRID '-' thisSID '.mat'],'Ripples')
        disp([thisRID '-' thisSID ' is finished!'])
    end
end