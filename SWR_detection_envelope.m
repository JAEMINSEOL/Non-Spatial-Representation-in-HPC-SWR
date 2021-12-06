%% This program is for SWR detection by
% - ver1 written on 2018-5-1


function  EEG_lowenvelope = SWR_detection_envelope(MotherROOT, sessionID,envelope_parameter,cscMat,cscRUN)

% disp('--------SWR_detection_envelope--------');

    
 MotherROOT = 'D:\SWR Analysis';

DataROOT = [ MotherROOT '\SWR envelope'];
%% set variables

% define loadCSC variables
CSCfileTag = 'RateReduced_150-250filtered';
exportMODE = 0;                                                      % 1 for extraction of all features to export to ncs files / 0 for extraction of eeg and timestamps only
behExtraction = 1;                                                     % 1 for extraction of behavior epoch only / 0 for extraction of whole session

%%

        findhypen = strfind(sessionID, '-');
        thisRID = sessionID(1:findhypen(1)-1);
        thisSID = sessionID(findhypen(1)+1:end);
       
%       if exist([MotherROOT '\rat' thisRID '\rat' thisRID '-' thisSID '\CSC' num2str(cscRUN) '_RateReduced.ncs'])
        cscID = cscMat{cscRUN};
    
        %% load csc data
        CSCdata = loadCSC(cscID, MotherROOT, CSCfileTag, exportMODE, behExtraction);
        load ([DataROOT '\lowpass\lowpass_20hz_' thisRID '-' thisSID '-' num2str(cscRUN) '.mat']);

        %% reshape eeg and timestamps for trial parsing
        eeg_expand = EEG.lowfil(:);
        
        Timestamps_expand = CSCdata.Timestamps ./ 1000000; % make to sec unit
        clear CSCdata.eeg CSCdata.Timestamps
        
        for i = 1 : length(Timestamps_expand)
            for count = 2 : 512
                Timestamps_expand(count, i) = Timestamps_expand(1, i) + (1/envelope_parameter.Fs)*(count-1);
            end
        end
        Timestamps_expand = Timestamps_expand(:);
        
    %% get Gaussian-smoothed envelope of eeg
        
        envelope = abs(hilbert(eeg_expand));
        envelope_smoothed = smoothdata(envelope, 'gaussian', envelope_parameter.gaussianSTD);
        
    %% remove noise period using session threshold
        temp_stat(1) = mean(envelope_smoothed);
        temp_stat(2) = std(envelope_smoothed,1);
        
        envelope_parameter.threshold(1) = temp_stat(1)+ envelope_parameter.noiseSTD * temp_stat(2);
        
        aboveNoiseThreshold = find(envelope_smoothed > envelope_parameter.threshold(1));
        
        envelope_smoothed_noiseRemoved = envelope_smoothed;
        envelope_smoothed_noiseRemoved(aboveNoiseThreshold,1) = NaN;
    
    %
        envelope_parameter.envelope_stat(1) = nanmean(envelope_smoothed_noiseRemoved);
        envelope_parameter.envelope_stat(2) = nanstd(envelope_smoothed_noiseRemoved,1);
        
        envelope_parameter.threshold(2) = envelope_parameter.envelope_stat(1) + ...
            envelope_parameter.thresholdSTD * envelope_parameter.envelope_stat(2);
    
    %% ripple detection
        ripples = []; ripples_index = []; % ripples_index is required for 3TToverlap

        aboveThreshold = find((envelope_smoothed<envelope_parameter.threshold(1))&(envelope_smoothed > envelope_parameter.threshold(2)));

    %%
        
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
            ripples(:,3) = ripples(:,2) - ripples(:,1);
        
        % ripples more than minimum duration
            ripples_index(ripples(:,3) < envelope_parameter.minDuration, :) = [];
            ripples(ripples(:,3) < envelope_parameter.minDuration, :) = [];
        
        % grouping
            temp = []; % contain row indices to be deleted
            for i = 2 : size(ripples,1)
                if ripples(i,1) - ripples(i-1,2) <= envelope_parameter.groupingInterval
                    ripples(i-1,2) = ripples(i,2);
                    ripples_index(i-1,2) = ripples_index(i,2);
                    temp(end+1) = i;
                end
            end
            ripples_index(temp,:) = [];
            ripples(temp,:) = [];
            
            ripples(:,3) = [];
            
            % remove noise ripples
            for i = 1 : size(ripples,1)
                if and(envelope_smoothed(ripples_index(i,1))>envelope_smoothed(ripples_index(i,1)-1), envelope_smoothed(ripples_index(i,2))<envelope_smoothed(ripples_index(i,2)+1))
                    ripples_index(i,1) = NaN; ripples_index(i,2) = NaN;
                    ripples(i,1) = NaN; ripples(i,2) = NaN;
                elseif and(envelope_smoothed(ripples_index(i,1))<envelope_smoothed(ripples_index(i,1)-1), envelope_smoothed(ripples_index(i,2))>envelope_smoothed(ripples_index(i,2)+1))
                    ripples_index(i,1) = NaN; ripples_index(i,2) = NaN;
                    ripples(i,1) = NaN; ripples(i,2) = NaN;
                end
            end
            
            ripples_index(isnan(ripples_index(:,:)))=[];
            ripples(isnan(ripples(:,:)))=[];
            if size(ripples_index,2) ~= 2
                ripples_index=reshape(ripples_index, [length(ripples_index)/2,2]);
                ripples=reshape(ripples, [length(ripples)/2,2]);
            end
        
        else
        disp('aboveThreshold is empty');
        end
    
    %% save as mat file
    
%         if ~exist([DataROOT '\rat' thisRID '\rat' thisRID '-' thisSID '\' cscID '.mat'], 'file')
%             save([DataROOT '\rat' thisRID '\rat' thisRID '-' thisSID '\' cscID '.mat'], 'envelope_parameter', 'ripples', 'ripples_index');
%         elseif exist([DataROOT '\rat' thisRID '\rat' thisRID '-' thisSID '\' cscID '.mat'], 'file')
%             save([DataROOT '\rat' thisRID '\rat' thisRID '-' thisSID '\' cscID '.mat'], 'envelope_parameter', 'ripples', 'ripples_index', '-append');
%         end
        EEG_lowenvelope.envelope_parameter = envelope_parameter;
        EEG_lowenvelope.envelope = envelope;
        EEG_lowenvelope.envelope_smoothed = envelope_smoothed;
%         ratst.envelope_smoothed_noiseRemoved= envelope_smoothed_noiseRemoved;
        EEG_lowenvelope.ripples = ripples;
        EEG_lowenvelope.ripples_index = ripples_index;
        
%          SWRenvelope.(['rat' thisRID '_' thisSID '_' num2str(cscRUN)])= ratst;
    %%
         clear CSCdata eeg_expand Timestamps_expand envelope envelope_smoothed aboveThreshold ripples
%     
        disp([cscID ' SWR detection envelope : done']);
    end
    
%       % cscRUN
% 
% 
% 
%   save('D:\SWR Analysis\SWRenvelope_1006_2std.mat','SWRenvelope')

