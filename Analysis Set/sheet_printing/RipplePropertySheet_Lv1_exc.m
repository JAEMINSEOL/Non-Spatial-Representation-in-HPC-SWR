function RipplePropertySheet_Lv1_exc(ROOT,SList,ripples,EEG,TargetTT,Spike,clusters,Params_Ripple)
unit = 5.0000e-04;
EEGpos = [.05 .95 .9 .03]; 
microSEC = 1e-06;
len = 10000;


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

epochST(2) = (epochST(1) - Tori) / unit;
epochED(2) = (epochED(1) - Tori) / unit;

nTT = size(TargetTT,1);

thisEPOCH(:,1) = [epochST(2):len:epochED(2)];
thisEPOCH(:,2) = [epochST(1):unit*len:epochED(1)];
%%
ti = 2;
for epoc = 1: size(thisEPOCH,1)-1
    sheet = figure;
    subplot('position', [.0 .96 .95 .03]);
    D = datetime;
    text(0.02, 0.5, ['Rat' RID '-' SID '-' SList.experimenter{1} '-' SList.type{1} '-d' num2str(SList.day)],'FontSize',14);
    text(0.85, 0.5, ['printed on ' datestr(D,'mmm-dd-yyyy')],'FontSize',10);
    axis off
    
    sheetPOS = [2 2 29.7 21.0];
    set(gcf,'units','centimeters','position',sheetPOS,'Color', [1 1 1]);
    thisRIPs = find(ripples.Var3>thisEPOCH(epoc,2) & ripples.Var3<=thisEPOCH(epoc+1,2));
    if isempty(thisRIPs)
        close all
        continue;
    end
    
    for t=1:nTT
        try
            y = EEG.(['TT' num2str(TargetTT(t))]).Raw(thisEPOCH(epoc,1):thisEPOCH(epoc+1,1));
            subplot('position', EEGpos - [0 .03 0 0]*t)
            plot(y,'k')
            xlim([0 size(y,1)])
            axis off
            if ~isempty(thisRIPs)
                for rip = 1:size(thisRIPs,1)
                    if ripples.ensemble(thisRIPs(rip))>=3, c='r'; else, c='k'; end
                    x = [ripples.Var1(thisRIPs(rip)), ripples.Var2(thisRIPs(rip))] - thisEPOCH(epoc,1);
                    p = patch([x(1) x(2) x(2) x(1)], [min(y) min(y) max(y) max(y)],c);
                    p.FaceAlpha = 0.3;
                    p.EdgeAlpha = 0;
                end
            end
        end
    end
            %% remove noise period using session threshold
            envelope_smoothed = EEG.(['TT' num2str(TargetTT(t))]).Envelope_smoothed(epochST(2):epochED(2));
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
            y = EEG.(['TT' num2str(TargetTT(t))]).Filtered(thisEPOCH(epoc,1):thisEPOCH(epoc+1,1));
            subplot('position', EEGpos - [0 .03 0 0]*t)
            plot(y,'b')
            line([0 len],[Params_Ripple.threshold(2) Params_Ripple.threshold(2)],'color','k')
            xlim([0 size(y,1)])
            axis off
            if ~isempty(thisRIPs)
                for rip = 1:size(thisRIPs,1)
                    if ripples.ensemble(thisRIPs(rip))>=3, c='r'; else, c='k'; end
                    x = [ripples.Var1(thisRIPs(rip)), ripples.Var2(thisRIPs(rip))] - thisEPOCH(epoc,1);
                    p = patch([x(1) x(2) x(2) x(1)], [min(y) min(y) max(y) max(y)],c);
                    p.FaceAlpha = 0.3;
                    p.EdgeAlpha = 0;
                end
            end
    
    %%
    Spkpos = EEGpos - [0 .03 0 0]*(t+1);

    
    spks_all=[];u=0; Units={};
    for un = 1:size(clusters,1)
        thisTTID = num2str(clusters.TT(un));
        thisCLID = num2str(str2double(clusters.ID{un}(end-1:end)));
        Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);

        thisSpks = Spk.t_spk(Spk.t_spk>thisEPOCH(epoc,2) & Spk.t_spk<=thisEPOCH(epoc+1,2));
        if ~isempty(thisSpks)
            spks_all = [spks_all;[thisSpks,ones(size(thisSpks,1),1)*u]];
            u=u+1;
            Units = [Units; [thisTTID '-' thisCLID]];
        end
    end
    subplot('position',Spkpos -[0 .17 0 -.12])
    hold on
    
    for s=1:size(spks_all,1)
        x = (spks_all(s,1) - thisEPOCH(epoc,2))/unit;
        patch([x x+ti x+ti x], [spks_all(s,2)+.2 spks_all(s,2)+.2 spks_all(s,2)+.8 spks_all(s,2)+.8],'k')
    end
    for rip = 1:size(thisRIPs,1)
        if ripples.ensemble(thisRIPs(rip))>=3, c='r'; else, c='k'; end
        x = [ripples.Var1(thisRIPs(rip)), ripples.Var2(thisRIPs(rip))] - thisEPOCH(epoc,1);
        p = patch([x(1) x(2) x(2) x(1)], [0 0 u u],c);
        p.FaceAlpha = 0.3;
        p.EdgeAlpha = 0;
    end
    xlim([0 size(y,1)])
    xticks([0:1000:len])
    xticklabels([thisEPOCH(epoc,2)-epochST:1000*unit:thisEPOCH(epoc+1,2)-epochST])
    xlabel('time(s)')
    
    yticks([0.5:1:u-0.5])
    yticklabels(Units)
    set(gca,'fontsize',5)
    
   cd(ROOT.Save)
saveImage(sheet, ['rat' RID '-' SID '_' num2str(epoc) '.png'], 'centimeters', sheetPOS); 
end


end
