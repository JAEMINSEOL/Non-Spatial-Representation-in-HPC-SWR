Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R2'];
ROOT.Fig1 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R1'];
ROOT.Fig2 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R2'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U0'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
if ~exist(ROOT.Fig1), mkdir(ROOT.Fig1); end
if ~exist(ROOT.Fig2), mkdir(ROOT.Fig2); end

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
Cluster_List = readtable([ROOT.Units '\ClusterList_SWR_CA1_filtered.xlsx']);


thisRegion = 'CA1';
RipplesTable_ori = readtable([ROOT.Old '\RipplesTable_Behav_' thisRegion '.xlsx']);
RipplesTable = readtable([ROOT.Old '\RipplesTable_Behav_' thisRegion '_speed_filtered.xlsx'],'ReadRowNames',false);
RipplesTable = RipplesTable_ori(RipplesTable_ori.ensemble>=3,:);

Experimenter = {'LSM','JS','SEB'};
ripple = [knnsearch(f,Params.Ripple(1)) knnsearch(f,Params.Ripple(2))];
unit = 5.0000e-04; ti = 0.5;

fd = dir(ROOT.Old);

RipplesTable_all = table;
thisSID_old = '';


for sid=1:size(RipplesTable,1)
    thisRip = RipplesTable(sid,:);
    thisSID = [jmnum2str(thisRip.rat,3) '-' jmnum2str(thisRip.session,2)];
    
    if ~strcmp(thisSID, thisSID_old)
        Recording_region_TT = Recording_region(thisSID,:);
        TargetTT = find(cellfun(cellfind(thisRegion),table2array(Recording_region_TT)'));
        EEG = LoadEEGData(ROOT, thisSID, TargetTT,Params,Params_Ripple);
        NumTT = length(fieldnames(EEG));
        Pos = load([ROOT.Raw.Mother '\rat' thisSID(1:3) '\rat' thisSID '\ParsedPosition.mat']);
        clusters = Cluster_List(Cluster_List.rat==thisRip.rat & Cluster_List.session==thisRip.session,:);
        Spike=LoadSpikeData(ROOT, thisSID, unique(clusters.TT), cellfindn);
        
        disp([thisSID ' plotting...'])
        thisSID_old = thisSID;
    end
    
    figure;
    % title
    subplot(4+fix(NumTT/6),6,1)
    title(thisRip.ID,'fontsize',15)
    axis off
    
    subplot(4+fix(NumTT/6),6,7)
    title(['trial ' (thisRip.trial{1}(end-2:end)) ', context ' num2str(thisRip.context)],'fontsize',15)
    axis off
    
    % position
    subplot(4+fix(NumTT/6),6,[5 6 17 18])
    if thisRip.speed>5 || thisRip.area~=0, cl='k'; else, cl='r'; end
    scatter(Pos.x,Pos.y,5,[.7 .7 .7],'filled')
    hold on
    scatter(thisRip.PosX,thisRip.PosY,20,cl,'filled')
    text(thisRip.PosX+3,thisRip.PosY,[jjnum2str(thisRip.speed,2) 'cm/s'],'color',cl)
    axis off
    
    % LFP trace
    t2=1; RipplePower=[];
    for t1=1:24
        subplot(4+fix(NumTT/6),6,18+t2)
        
        try
            axis off
            plot(EEG.(['TT' num2str(t1)]).Raw(thisRip.STindex:thisRip.EDindex))
            hold on
            plot(EEG.(['TT' num2str(t1)]).Filtered(thisRip.STindex:thisRip.EDindex))
            title(['TT' num2str(t1)])
            [pxx, f] = CalPSD(EEG.(['TT' num2str(t1)]).Raw(thisRip.STindex:thisRip.EDindex), 'k',[0:300],0.01,0);
    RipplePower(t1) = nanmean(10*log10(pxx(ripple(1):ripple(2))));
    RipplePower(t1) = nanmean(EEG.(['TT' num2str(t1)]).Filtered(thisRip.STindex:thisRip.EDindex));
            axis off
            t2 = t2+1;
        end
    end

    
        % reactivated cells
    subplot(4+fix(NumTT/6),6,[3,4])
    [~,m] = max(RipplePower);
     thisEEG = EEG.(['TT' num2str(m)]).Raw(thisRip.STindex:thisRip.EDindex);
            plot(thisEEG)
            hold on
            plot(EEG.(['TT' num2str(m)]).Filtered(thisRip.STindex:thisRip.EDindex))
            title(['TT' num2str(m)])
                x1=0; x2=thisRip.RippleDuration*2e3;
    xlim([x1-20 x2+20])
    line([x1 x1], [min(thisEEG) max(thisEEG)+50], 'color','k','linestyle','--')
    line([x2 x2], [min(thisEEG) max(thisEEG)+50], 'color','k','linestyle','--')
            axis off
            
            
    subplot(4+fix(NumTT/6),6,[9,10,15,16])
    
    spks_epoch=[];u=0; Units={}; 
    cls = size(clusters,1);
    for un = 1:cls
        thisTTID = num2str(clusters.TT(un));
        thisCLID = num2str(str2double(clusters.ID{un}(end-1:end)));
        Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);
        
        thisSpks = Spk.t_spk(Spk.t_spk>thisRip.STtime & Spk.t_spk<=thisRip.EDtime);
        if ~isempty(thisSpks)
            spks_epoch = [spks_epoch;[thisSpks,ones(size(thisSpks,1),1)*un,ones(size(thisSpks,1),1)*u]];
            u=u+1;
            
        end
        Units = [Units; [thisTTID '-' thisCLID]];
    end

    hold on
    for s=1:size(spks_epoch,1)
        x = (spks_epoch(s,1) - thisRip.STtime)*1e3;
        patch([x-ti x+ti x+ti x-ti], [spks_epoch(s,3)+.2 spks_epoch(s,3)+.2 spks_epoch(s,3)+.8 spks_epoch(s,3)+.8],'k')
    end
    x1=0; x2=thisRip.RippleDuration*1e3;
    xlim([x1-10 x2+10])
    line([x1 x1], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
    line([x2 x2], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
    axis off

end