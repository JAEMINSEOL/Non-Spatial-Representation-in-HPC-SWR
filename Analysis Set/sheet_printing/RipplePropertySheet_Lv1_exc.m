function RipplePropertySheet_Lv1_exc(ROOT,ID,ripples,EEG,TargetTT,Spike,clusters)
unit = 5.0000e-04;
EEGpos = [.05 .95 .9 .03]; 
microSEC = 1e-06;
len = 10000;

RID = ID(1:3);
SID = ID(5:6);

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
    text(0.02, 0.5, ['Rat' RID '-S' SID ],'FontSize',14);
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
        subplot('position', EEGpos - [0 .03 0 0]*t)
        y = EEG.(['TT' num2str(TargetTT(t))]).Raw(thisEPOCH(epoc,1):thisEPOCH(epoc+1,1));
        plot(y,'k')
        xlim([0 size(y,1)])
        axis off
        if ~isempty(thisRIPs)
            for rip = 1:size(thisRIPs,1)
                x = [ripples.Var1(thisRIPs(rip)), ripples.Var2(thisRIPs(rip))] - thisEPOCH(epoc,1);
                p = patch([x(1) x(2) x(2) x(1)], [min(y) min(y) max(y) max(y)],'r');
                p.FaceAlpha = 0.3;
                p.EdgeAlpha = 0;
            end
        end
    end

    Spkpos = EEGpos - [0 .03 0 0]*t;

    
    spks_all=[];u=0; Units={};
    for un = 1:size(clusters,1)
        thisTTID = num2str(clusters.TT(un));
        thisCLID = num2str(str2double(clusters.ID{un}(end-1:end)));
        Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);
        thisSpks = Spk.t_spk(Spk.t_spk>thisEPOCH(epoc,2) & Spk.t_spk<=thisEPOCH(epoc+1,2));
        if ~isempty(thisSpks)
            spks_all = [spks_all;[thisSpks,ones(size(thisSpks,1),1)*u]];
            u=u+1;
            Units = [Units; ['TT' thisTTID '-' thisCLID]];
        end
    end
    subplot('position',Spkpos -[0 .2 0 -.15])
    hold on
    for s=1:size(spks_all,1)
        x = (spks_all(s,1) - thisEPOCH(epoc,2))/unit;
        patch([x x+ti x+ti x], [spks_all(s,2)+.2 spks_all(s,2)+.2 spks_all(s,2)+.8 spks_all(s,2)+.8],'k')
    end
    for rip = 1:size(thisRIPs,1)
        x = [ripples.Var1(thisRIPs(rip)), ripples.Var2(thisRIPs(rip))] - thisEPOCH(epoc,1);
        p = patch([x(1) x(2) x(2) x(1)], [0 0 u u],'r');
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
saveImage(sheet, ['rat' ID '_' num2str(epoc) '.png'], 'centimeters', sheetPOS); 
end


end
