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
            
            %% ripple detection
            ripples = []; ripples_index = []; % ripples_index is required for 3TToverlap
            
            aboveThreshold = find((envelope_smoothed<Params_Ripple.threshold(1))&(envelope_smoothed > Params_Ripple.threshold(2)));
            
            if ~isempty(aboveThreshold)
                ripples_index(end+1, 1) = aboveThreshold(1,1);
                for i = 2 : length(aboveThreshold)
                    if aboveThreshold(i,1) - aboveThreshold(i-1,1) > 1
                        ripples_index(end, 2) = aboveThreshold(i-1,1);
                        ripples_index(end+1, 1) = aboveThreshold(i,1);
                    end
                end
                ripples_index(end,2) = aboveThreshold(end,1);
                
                % index to timestamp
                ripples(:,1) = Timestamps_expand(ripples_index(:,1));
                ripples(:,2) = Timestamps_expand(ripples_index(:,2));

                
                % grouping
                temp = []; % contain row indices to be deleted
                ripples2=[]; ripples_index2=[]; 
                i=1;
                while i<size(ripples,1)
                    j=i+1;
                        while ripples(j,1) - ripples(j-1,2) <= Params_Ripple.groupingInterval
                            j=j+1;
                            if j>size(ripples,1), continue; end
                        end
                    ripples2 = [ripples2;[ripples(i,1),ripples(j-1,2)]];
                    ripples_index2 = [ripples_index2;[ripples_index(i,1),ripples_index(j-1,2)]];
                    i=j;
                end
                if ripples_index(end,2)~=ripples_index2(end,2)
                    ripples_index2 = [ripples_index2 ; ripples_index(end,1:2)];
                    ripples2 = [ripples2 ; ripples(end,1:2)];
                end
                ripples_index= ripples_index2;
                ripples = ripples2;
                
                                
                
                % ripples more than minimum duration
                ripples(:,3) = ripples(:,2) - ripples(:,1);
                ripples_index(ripples(:,3) < Params_Ripple.minDuration, :) = [];
                ripples(ripples(:,3) < Params_Ripple.minDuration, :) = [];
                ripples(:,3)=[];
                % remove noise ripples
%                 for i = 1 : size(ripples,1)
%                     if and(envelope_smoothed(ripples_index(i,1))>envelope_smoothed(ripples_index(i,1)-1), envelope_smoothed(ripples_index(i,2))<envelope_smoothed(ripples_index(i,2)+1))
%                         ripples_index(i,1) = NaN; ripples_index(i,2) = NaN;
%                         ripples(i,1) = NaN; ripples(i,2) = NaN;
%                     elseif and(envelope_smoothed(ripples_index(i,1))<envelope_smoothed(ripples_index(i,1)-1), envelope_smoothed(ripples_index(i,2))>envelope_smoothed(ripples_index(i,2)+1))
%                         ripples_index(i,1) = NaN; ripples_index(i,2) = NaN;
%                         ripples(i,1) = NaN; ripples(i,2) = NaN;
%                     end
%                 end
                
                ripples_index(isnan(ripples_index(:,:)))=[];
                ripples(isnan(ripples(:,:)))=[];

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
