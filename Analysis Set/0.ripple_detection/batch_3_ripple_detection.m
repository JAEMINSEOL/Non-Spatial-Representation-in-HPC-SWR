Initial_SWRFilter_common;
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

Experimenter = {'LSM','JS','SEB'};



fd = dir(ROOT.Save);

RipplesTable_all = table;

for sid=1:size(SessionList,1)
    if SessionList.include(sid) && ismember(SessionList.experimenter(sid),Experimenter)
        thisRID = jmnum2str(SessionList.rat(sid),3);
        thisSID = jmnum2str(SessionList.session(sid),2);
        
        TargetTT = [1:24]';
        
        %% EEG
        EEG = LoadEEGData(ROOT, [thisRID '-' thisSID], TargetTT,Params,Params_Ripple);
        Ripples = struct;
        for thisTTID=1:24
            try
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
            
            aboveThreshold = (envelope_smoothed<Params_Ripple.threshold(1))&(envelope_smoothed > Params_Ripple.threshold(2));
            beginThreshold =(envelope_smoothed > Params_Ripple.threshold(3));
            
            
            if ~isempty(aboveThreshold)
                ripples_index = zeros(length(aboveThreshold),1);
                
                for i=1:length(aboveThreshold)-1
                        if aboveThreshold(i)==0 
                            if aboveThreshold(i+1)==1
                                j=i;
                                while beginThreshold(j)==1
                                    aboveThreshold(j)=1;
                                    ripples_index(j)=1;
                                    j=j-1;
                                    if j<1, break; end
                                end
                            elseif aboveThreshold(i+1)==0
                                ripples_index(i)=0;
                            end
                        elseif aboveThreshold(i)==1
                            if aboveThreshold(i+1)==0
                                j=i;
                                while beginThreshold(j)==1
                                    aboveThreshold(j)=1;
                                    ripples_index(j)=1;
                                    j=j+1;
                                    if j>length(aboveThreshold), break; end
                                end
                            elseif aboveThreshold(i+1)==1
                                ripples_index(i)=1;
                            end
                        end
                end
                                              
                ripples = zeros(1,1);
                ripples(:,1) = Timestamps_expand(end);
                ripples(:,2) = EEG.(['TT' num2str(thisTTID)]).Raw(end);
                ripples(:,3) = envelope_smoothed_noiseRemoved(end);

                Ripples.(['TT' num2str(thisTTID)]).ripples = ripples;
                Ripples.(['TT' num2str(thisTTID)]).index = ripples_index;
                
                disp([thisRID '-' thisSID '-' num2str(thisTTID) ', is finished!']);
            else
                disp([thisRID '-' thisSID '-' num2str(thisTTID) ', aboveThreshold is empty']);
            end
            catch
            end
        end
        save([ROOT.Save '\' thisRID '-' thisSID '.mat'],'Ripples')
        disp([thisRID '-' thisSID ' is finished!'])
    end
end
