function RipplePropertySheet_Lv1_exc(ROOT,SList,ripples,EEG,TargetTT,Spike,clusters,Params_Ripple)
unit = 5.0000e-04;
EEGpos = [.05 .96 .9 .025];
microSEC = 1e-06;
len = 20000;


RID = jmnum2str(SList.rat,3);
SID = jmnum2str(SList.session,2);
%Load Epoch information
cd([ROOT.Raw.Mother '\rat' RID]);

if exist(['behaviorEpoch_rat' RID '.xlsx'])
    EPOCH = xlsread(['behaviorEpoch_rat' RID '.xlsx']);
else
    EPOCH = csvread(['behaviorEpoch_rat' RID '.csv']);
end
Tori = ripples.Var3(1) - ripples.Var1(1) * unit;

epochST(1) = EPOCH(str2double(SID),1) * microSEC;
epochED(1) = EPOCH(str2double(SID),2) * microSEC;

epochST(2) = floor((epochST(1) - Tori) / unit+1);
epochED(2) = ceil((epochED(1) - Tori) / unit);

nTT = size(TargetTT,1);

thisEPOCH(:,1) = [epochST(2):len:epochED(2)];
thisEPOCH(:,2) = [epochST(1):unit*len:epochED(1)];
%%
ti = 2;
filt_eeg_all=[];
envelope_smoothed=[];
for t=1:nTT
    try
        filt_eeg_all = [filt_eeg_all,EEG.(['TT' num2str(TargetTT(t))]).Filtered];
        envelope_smoothed = [envelope_smoothed,EEG.(['TT' num2str(TargetTT(t))]).Envelope_smoothed];
    end
    if t==nTT
        envelope_smoothed = mean(envelope_smoothed,2);
        filt_eeg_all = mean(filt_eeg_all,2);
        temp_stat(1) = mean(envelope_smoothed);
        temp_stat(2) = std(envelope_smoothed,1);
        
        Params_Ripple.threshold(1) = temp_stat(1)+ Params_Ripple.noiseSTD * temp_stat(2);
        aboveNoiseThreshold = find(envelope_smoothed > Params_Ripple.threshold(1));
        envelope_smoothed_noiseRemoved = envelope_smoothed;
        envelope_smoothed_noiseRemoved(aboveNoiseThreshold,1) = NaN;

        Params_Ripple.envelope_stat(1) = nanmean(envelope_smoothed_noiseRemoved);
        Params_Ripple.envelope_stat(2) = nanstd(envelope_smoothed_noiseRemoved,1);
        
        Params_Ripple.threshold(2) = Params_Ripple.envelope_stat(1) + Params_Ripple.thresholdSTD * Params_Ripple.envelope_stat(2);
    end
end
for epoc = 1: size(thisEPOCH,1)-1
    sheet = figure;
    subplot('position', [.0 .97 .95 .02]);
    D = datetime;
    text(0.02, 0.5, ['Rat' RID '-' SID '-' SList.experimenter{1} '-' SList.type{1} '-d' num2str(SList.day) '-epoch' num2str(epoc) '-60Hz filtered'],'FontSize',14);
    text(0.85, 0.5, ['printed on ' datestr(D,'mmm-dd-yyyy')],'FontSize',10);
    axis off
    
    sheetPOS = [20 20 2000 1500];
    set(gcf,'position',sheetPOS,'Color', [1 1 1]);
    thisRIPs = find((ripples.Var3>thisEPOCH(epoc,2) & ripples.Var3<=thisEPOCH(epoc+1,2)) |...
        (ripples.Var4>thisEPOCH(epoc,2) & ripples.Var4<=thisEPOCH(epoc+1,2)));
    if isempty(thisRIPs)
        close all
        continue;
    end

    for t=1:nTT
        try
            
            y = EEG.(['TT' num2str(TargetTT(t))]).Raw(thisEPOCH(epoc,1):thisEPOCH(epoc+1,1));
            subplot('position', EEGpos - [0 .025 0 0]*t)
            
            % 60Hz notch
            Fn = Params_Ripple.Fs/2;
            n=5;
            Wn = [55 65];
            ftype = 'stop';
            [b,a] = butter(n,Wn/Fn,ftype);
            y = filter(b,a,y);
            
            plot(y,'k')
            xlim([0 size(y,1)])
            axis off
            if ~isempty(thisRIPs)
                for rip = 1:size(thisRIPs,1)
                    if ripples.Var5(thisRIPs(rip))>=3
                        c='r';
                    else, c='k'; end
                    x = [ripples.Var1(thisRIPs(rip)), ripples.Var2(thisRIPs(rip))] - thisEPOCH(epoc,1);
                    p = patch([x(1) x(2) x(2) x(1)], [min(y) min(y) max(y) max(y)],c);
                    p.FaceAlpha = 0.3;
                    p.EdgeAlpha = 0;
                end
            end
            
            if t==nTT
                filt_eeg = filt_eeg_all(thisEPOCH(epoc,1):thisEPOCH(epoc+1,1));
                subplot('position', EEGpos - [0 .025 0 0]*t)
                plot(filt_eeg,'b')
                line([0 len],[Params_Ripple.threshold(2) Params_Ripple.threshold(2)],'color','k')
                xlim([0 size(filt_eeg ,1)])
                axis off
                if ~isempty(thisRIPs)
                    for rip = 1:size(thisRIPs,1)
                        if ripples.Var5(thisRIPs(rip))>=3, c='r'; else, c='k'; end
                        x = [ripples.Var1(thisRIPs(rip)), ripples.Var2(thisRIPs(rip))] - thisEPOCH(epoc,1);
                        p = patch([x(1) x(2) x(2) x(1)], [min(filt_eeg ) min(filt_eeg ) max(filt_eeg ) max(filt_eeg )],c);
                        p.FaceAlpha = 0.3;
                        p.EdgeAlpha = 0;
                    end
                end
            end
        end
    end
    
    
    %%
    Spkpos = EEGpos - [0 .025 0 0]*(t+1);
    
    
    spks_epoch=[];u=0; Units={};
    cls = size(clusters,1);
    for un = 1:cls
        thisTTID = num2str(clusters.TT(un));
        thisCLID = num2str(str2double(clusters.ID{un}(end-1:end)));
        Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);
        
        thisSpks = Spk.t_spk(Spk.t_spk>thisEPOCH(epoc,2) & Spk.t_spk<=thisEPOCH(epoc+1,2));
        if ~isempty(thisSpks)
            spks_epoch = [spks_epoch;[thisSpks,ones(size(thisSpks,1),1)*un]];
            
            
        end
        Units = [Units; [thisTTID '-' thisCLID]];
    end
    subplot('position',Spkpos -[0 .33 0 -.3])
    hold on
    if isempty(spks_epoch), close all; continue; end
    for s=1:size(spks_epoch,1)
        x = (spks_epoch(s,1) - thisEPOCH(epoc,2))/unit;
        patch([x x+ti x+ti x], [spks_epoch(s,2)+.2 spks_epoch(s,2)+.2 spks_epoch(s,2)+.8 spks_epoch(s,2)+.8],'k')
    end
    for rip = 1:size(thisRIPs,1)
        if ripples.Var5(thisRIPs(rip))>=3, c='r'; else, c='k'; end
        x = [ripples.Var1(thisRIPs(rip)), ripples.Var2(thisRIPs(rip))] - thisEPOCH(epoc,1);
        p = patch([x(1) x(2) x(2) x(1)], [0 0 cls+1 cls+1],c);
        p.FaceAlpha = 0.3;
        p.EdgeAlpha = 0;
        text(x(1),cls+2,num2str(thisRIPs(rip)),'color',c)
    end
    xlim([0 size(y,1)])
    xticks([0:1000:len])
    xticklabels([thisEPOCH(epoc,2)-epochST:1000*unit:thisEPOCH(epoc+1,2)-epochST])
    xlabel('time(s)')
    
    yticks([1.5:1:cls+0.5])
    yticklabels(Units)
    ylim([1, cls+1])
    set(gca,'fontsize',8)
 
    %%
    cd(ROOT.Fig)
    if ~exist(SList.experimenter{1})
        mkdir(SList.experimenter{1});
    end
    cd(SList.experimenter{1})
    if ~exist(RID)
        mkdir(RID);
    end
    cd(RID)
    saveas(sheet, ['rat' RID '-' SID '_' num2str(epoc) '.png']);
    close all
    clear sheet
end

clear all
end
