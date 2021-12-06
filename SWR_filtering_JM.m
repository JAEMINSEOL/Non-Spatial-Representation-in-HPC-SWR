
clear all; clc; fclose all;

ROOT.Mother = 'D:\HPC-SWR project';
ROOT.Raw.Mother = 'G:\EPhysRawData\RawData';
ROOT.Program = [ROOT.Mother '\Analysis Program'];
ROOT.Info = [ROOT.Mother '\Information Sheet'];
ROOT.Save = [ROOT.Mother '\Processed Data\Ripples\S4N100G40 - not grouped'];
addpath(genpath(ROOT.Program))
%% set initial parameters
% define frequency
Params.Fs = 2000;  % sampling frequency
Params.Fn = Params.Fs/2;  % Nyquist frequency
Params.F0 = 60/Params.Fn; % notch frequency

% define loadCSC variable
Params.CSCfileTag = 'RateReduced';
Params.exportMODE = 1;
Params.behExtraction = 0;

% switch
Params.noiseFilteringSwitch = 1;
Params.saveSwitch = 1;
Params.exportSwitch = 1;

% define data duration
Params.freqN = 2048;
Params.freqLimit = 450;
Params.freqBin = ceil(Params.freqLimit/Params.Fs * Params.freqN);

Params.noise = [350 450]; %noise Range
Params.Ripple = [150 250]; %SWR range
Params.low = 20;




%% load data list
DataList = readtable([ROOT.Info '\SessionList_SWR.xlsx']);

%% for LOOP
for clRUN=146:length(DataList.rat)
    try
    if ~DataList.include(clRUN), continue; end
    
    % set local parameter
    thisRID.n = DataList.rat(clRUN); if thisRID.n<100, thisRID.s=['0' num2str(thisRID.n)]; else thisRID.s=num2str(thisRID.n); end
    thisSID.n = DataList.session(clRUN); if thisSID.n<10, thisSID.s=['0' num2str(thisSID.n)]; else thisSID.s=num2str(thisSID.n); end
    thisSID.full = [thisRID.s '-' thisSID.s];
    thisStype = DataList.type(clRUN);
%     Params.CSCfileTag = cell2mat(DataList.CSCfileTag(clRUN));
Params.CSCfileTag = 'RateReduced';
    ROOT.Raw.Run = [ROOT.Raw.Mother '\rat' thisRID.s '\rat' thisSID.full];
    TTAvailable{1,clRUN} = thisSID.full;
    
    for thisTTID = 1:24
        cscID = [thisSID.full '-' num2str(thisTTID)];
        
        LocalP = exist([ROOT.Raw.Run '\' 'CSC' num2str(thisTTID) '_' Params.CSCfileTag '.ncs']);
        TTAvailable{thisTTID+1,clRUN} = ~(~LocalP);
        if ~LocalP, continue; end
        
        % load CSC data
        cscData = loadCSC(cscID, ROOT.Raw.Mother, Params.CSCfileTag, Params.exportMODE, Params.behExtraction);
%         cscData = loadCSC('000-00-1', ROOT.Raw.Mother, Params.CSCfileTag, Params.exportMODE, Params.behExtraction);
        [eeg_expand,Timestamps_expand] = expandCSC(cscData);
        EEG.Raw = eeg_expand(:);
        CSCdata_filtered.eeg = zeros(size(EEG.Raw));
        
        % 150-250Hz filtering
        Params.cRange = Params.Ripple;
        Params_Ripple = SetRippleParams(Params);
        LocalP = [Params.CSCfileTag '_' num2str(Params.cRange(1)) '-' num2str(Params.cRange(2)) 'filtered'];
        EEG.Filtered = FiltLFP(EEG.Raw,Params.Ripple,Params.Fn,'bandpass');
        
        % export filtered LFP (optional)
%         if ~exist([ROOT.Raw.Run '\' 'CSC' num2str(thisTTID) '_' LocalP '.ncs'])
if 1
                        % reshape
                        colN = floor(numel(EEG.Filtered)/512);
                        cscData.eeg = reshape(EEG.Filtered, 512, colN);
                        findHYPHEN = find(cscID == '\');
            
                        % export
                        Switch = 'butterSwitch_export';
                        filename_tail = LocalP;
                        export_Mat2NlxCSC(cscData, cscData, ROOT.Raw.Run, filename_tail, Switch)
                        clear filename_tail Switch
        end
        %% get Gaussian-smoothed envelope of eeg
        
        EEG.Envelope = abs(hilbert(EEG.Filtered));
        EEG.Envelope_smoothed = smoothdata(EEG.Envelope, 'gaussian', Params_Ripple.gaussianSTD);
        
        %% ripple detection
        ripples = []; ripples_index = []; % ripples_index is required for 3TToverlap
        RippleStat.mean = nanmean(EEG.Envelope_smoothed);
        RippleStat.std = nanstd(EEG.Envelope_smoothed,1);
        Params_Ripple.SignalBeginT = RippleStat.mean + Params_Ripple.beginthresholdSTD * RippleStat.std;
        Params_Ripple.SignalT = RippleStat.mean + Params_Ripple.thresholdSTD * RippleStat.std;
        Params_Ripple.NoiseT = RippleStat.mean + Params_Ripple.noiseSTD * RippleStat.std;
        
        aboveThreshold = find((EEG.Envelope_smoothed>Params_Ripple.SignalBeginT));
        
        %%
        if ~isempty(aboveThreshold)
            ripples_index(1, 1) = aboveThreshold(1,1);
            j=1;
            for i = 2 : length(aboveThreshold)
                if aboveThreshold(i,1) - aboveThreshold(i-1,1) > 1
                    ripples_index(j, 2) = aboveThreshold(i-1,1);
                    ripples_index(j+1, 1) = aboveThreshold(i,1);
                    j=j+1;
                end
            end
            ripples_index(j,2) = aboveThreshold(i,1);
            
            
            % index to timestamp
            ripples(:,1) = Timestamps_expand(ripples_index(:,1));
            ripples(:,2) = Timestamps_expand(ripples_index(:,2));
            ripples(:,3) = ripples(:,2) - ripples(:,1);
            
            
            
            % grouping
            temp2 = []; % contain row indices to be deleted
            ripples_grouped=[];
            ripples_index_grouped=[];
            noise_index_grouped=[];
%             j=1; k=1; l=1;
%             for i = 1 : size(ripples,1)-1
%                 if ripples(i+1,1) - ripples(i,2) <= Params_Ripple.groupingInterval
%                     
%                                     ripples_index(i,2) = ripples_index(i+1,2);
%                     temp2(end+1) = i;
%                     j=j+1;
%                 else
%                     ripples_grouped(l,1) = ripples(k,1);
%                     ripples_grouped(l,2) = ripples(j,2);
%                     ripples_index_grouped(l,1) = ripples_index(k,1);
%                     ripples_index_grouped(l,2) = ripples_index(j,2);
%                     l=l+1;
%                     k=i+1;
%                     j=i+1;
%                 end
%                 
%             end
                                ripples_grouped= ripples;
                    ripples_index_grouped= ripples_index;
   
            
            
            % ripples more than minimum duration
            ripples_grouped(:,3) = ripples_grouped(:,2)-ripples_grouped(:,1);
            ripples_index_grouped(ripples_grouped(:,3) < Params_Ripple.minDuration, :) = [];
            ripples_grouped(ripples_grouped(:,3) < Params_Ripple.minDuration, :) = [];
%              ripples_grouped(ripples_grouped(:,3) > Params_Ripple.maxDuration, :) = [];
            ripples_grouped(:,3) = [];
            
            % remove noise ripples
            j=1;
            for i=1:size(ripples_grouped)
                if max(EEG.Envelope_smoothed(ripples_index_grouped(i,1):ripples_index_grouped(i,2))) > Params_Ripple.NoiseT ||...
                        max(EEG.Envelope_smoothed(ripples_index_grouped(i,1):ripples_index_grouped(i,2))) < Params_Ripple.SignalT
                    if  max(EEG.Envelope_smoothed(ripples_index_grouped(i,1):ripples_index_grouped(i,2))) > Params_Ripple.NoiseT
                        noise_index_grouped(j,:)=ripples_index_grouped(i,:);
                        j=j+1;
                    end
                    ripples_index_grouped(i,1)=NaN;
                end
            end
            ripples_grouped(find(isnan(ripples_index_grouped(:,1))),:)=[];
            ripples_index_grouped(find(isnan(ripples_index_grouped(:,1))),:)=[];
            
        else
            disp('aboveThreshold is empty');
        end
        Ripples.(['TT' num2str(thisTTID)]).ripples = ripples_grouped;
        Ripples.(['TT' num2str(thisTTID)]).index = ripples_index_grouped;
        Noise.(['TT' num2str(thisTTID)]).index = noise_index_grouped;
        
        %%
        %             s=Timestamps_expand-Timestamps_expand(1);
        %             n1=1001; n2=3001;
        %             figure;
        %
        %             subplot(3,1,1)
        %
        %             plot(EEG.Raw(n1:n2));
        %               title([thisSID.full ', ' num2str(s(n1)) 's ~' num2str(s(n2)) 's'])
        %             subplot(3,1,2)
        %             plot(EEG.Filtered(n1:n2));
        %             hold on
        %             plot(EEG.Envelope(n1:n2));
        %             subplot(3,1,3)
        %             plot(EEG.Envelope_smoothed(n1:n2));
        %             ylim([-100 100])
        %
        %%
        sz=64;
        
%         for j = 1:fix(size(ripples_grouped,1)/sz)+1
%             %             temp1=sortrows(randperm(size(ripples_grouped,1),100)');
%             temp1=linspace(sz*(j-1)+1,sz*j,sz);
% %             figure('position',[100,50,2200,1300]);
%             for i1=1:sz*2
%                 sz2=[100,100];
%                 
%                 if mod(i1,2)==1
%                     i2=fix(i1/2)+1;
%                     if size(ripples_index_grouped,1)<temp1(i2), break; end
%                     n1=ripples_index_grouped(temp1(i2),1); n2=ripples_index_grouped(temp1(i2),2);
%                     if n1-sz2(1)<1, sz2(1)=n1-1; end
%                     if n2+sz2(2)>length(EEG.Raw), sz2(2)=length(EEG.Raw)-n2; end
%                     temp2 = EEG.Raw(n1-sz2(1):n2+sz2(2)); c='k';
%                     
%                 else
%                     i2=fix(i1/2);
%                     n1=ripples_index_grouped(temp1(i2),1); n2=ripples_index_grouped(temp1(i2),2);
%                     if n1-sz2(1)<1, sz2(1)=n1-1; end
%                     if n2+sz2(2)>length(EEG.Raw), sz2(2)=length(EEG.Raw)-n2; end
%                     temp2 = EEG.Envelope_smoothed(n1-sz2(1):n2+sz2(2)); c='b';
%                 end
%                 
% %                 subplot(sqrt(sz),sqrt(sz)*2,i1)
%                 hold on
%                 if mod(i1,2)==1, title(temp1(i2)); text(sz2(1),min(temp2),[num2str((n2-n1)*0.5)],'color','r');
%                     if size(ripples_index_grouped,1)>temp1(i2), text(n2-n1+sz2(1),max(temp2)*0.9,[num2str((ripples_index_grouped(temp1(i2)+1,1)-n2)*0.5)],'color','m'); end
%                 else
% %                     line([0 n2-n1+sz2(1)+sz2(2)],[Params_Ripple.SignalT Params_Ripple.SignalT], 'color','m');
% %                         line([0 n2-n1+sz2(1)+sz2(2)],[Params_Ripple.SignalBeginT Params_Ripple.SignalBeginT], 'color','m');
%                         %                 line([0 ripples_index_grouped(i,2)-ripples_index_grouped(i,1)+200],[Params_Ripple.NoiseT Params_Ripple.NoiseT], 'color','k');
%                 end
%                 %                 plot(EEG.Filtered(ripples_index_grouped(i,1)-100:ripples_index_grouped(i,2)+100))
% %                 plot(temp2,'color',c)
% %                 line([sz2(1) sz2(1)], [min(temp2) max(temp2)], 'color','r')
% %                 line([n2-n1+sz2(1) n2-n1+sz2(1)], [min(temp2) max(temp2)], 'color','r')
%                 %                 line([0 ripples_index_grouped(i,2)-ripples_index_grouped(i,1)+200],[Params_Ripple.SignalT Params_Ripple.SignalT],'Color','red')
%                 %                 line([0 ripples_index_grouped(i,2)-ripples_index_grouped(i,1)+200],[Params_Ripple.NoiseT Params_Ripple.NoiseT],'Color','k')
%                 
%                 axis off
%                 
%             end
% %             sgtitle(['Signal>' num2str(Params_Ripple.thresholdSTD) '(onset>' num2str(Params_Ripple.beginthresholdSTD) ')std, Noise>' num2str(Params_Ripple.noiseSTD) ...
% %                 'std, Gaussian filtered(' num2str(Params_Ripple.gaussianSTD*0.5) 'ms), after 150-250 Hz bandpass filtering'])
% %             cd([ROOT.Save '\Plot'])
% %             if ~exist([ROOT.Save '\Plot\' thisSID.full]), mkdir(thisSID.full);  end
% %             saveas(gcf,[ROOT.Save '\Plot\' thisSID.full  '\' thisSID.full '-TT' num2str(thisTTID) '-' num2str(j) '.png'])
% % %             close gcf;
%         end
        disp([thisSID.full ' TT' num2str(thisTTID) ' is finished!'])
        
    end
    save([ROOT.Save '\' thisSID.full '.mat'], 'Ripples')
    save([ROOT.Save '\' thisSID.full '.mat'], 'Noise','-append')
    disp([thisSID.full ' has been done!']);
    catch
        disp([thisSID.full ' failed!']);
    end
end

%% ----


% clear CSCdata thiseeg filtered_eeg PSD_afterFiltering PSD_beforeFiltering f

function Params_Ripple = SetRippleParams(Params)
Params_Ripple.Fs = 2000;                                   % sampling frequency
Params_Ripple.gaussianSTD = 40;                            % gaussianSTD / Fs = moving window size for gaussian smoothing (sec) 0.5ms=1idx
Params_Ripple.thresholdSTD = 4;                            % SD
Params_Ripple.beginthresholdSTD = 1;                           % begin & end SD
Params_Ripple.noiseSTD = 100;                               % 11 SD
Params_Ripple.minDuration = 0.02;                          % sec
Params_Ripple.maxDuration = 0.1;    
Params_Ripple.groupingInterval = 0.02;                     % sec
Params_Ripple.boundarySTD = 1;                             % SD

if isequal(Params.cRange, Params.noise) %noise range
    Params_Ripple.thresholdSTD = 2;                            % SD
    Params_Ripple.noiseSTD = 11;                               %  not used for noise range
    Params_Ripple.groupingInterval = 2;                     % sec
end
end