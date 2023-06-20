Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Mother '\Processed Data\ripples_mat\R4'];
ROOT.Rip5 = [ROOT.Mother '\Processed Data\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Mother '\Manuscript figures\R0_fig'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U2'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];

dir = '';

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

unit = 5.0000e-04;
EEGpos = [.05 .96 .9 .035];
microSEC = 1e-06;
len = 20000;
thisFRMapSCALE=2;
%%
thisRegion = 'CA1';
thisRegion2 = 'CA1_field';
Sel = 'LScene';
RipplesTable = readtable([ROOT.Rip '\RipplesTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);

TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

for sid=64:size(SessionList,1)
    if ~SessionList.include(sid) | ~strcmp(SessionList.experimenter(sid),'LSM'), continue; end
    thisRIDn = SessionList.rat(sid); thisRID=jmnum2str(thisRIDn,3);
    thisSIDn = SessionList.session(sid); thisSID=jmnum2str(thisSIDn,2);
    thisRSID = [thisRID '-' thisSID];
    load([ROOT.Behav '\' thisRSID '.mat']);
    BehavTable = readtable([ROOT.Behav '\' thisRSID '.xlsx']);
    Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);


    Pos = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\ParsedPosition.mat']);
    diverging_point = get_divergingPoint(ROOT.Info, thisRID, thisSID);

    diverging_point = diverging_point*0.23;
    stem_end_index =  (max(Pos.y)-diverging_point)/thisFRMapSCALE;

    Recording_region_TT = Recording_region(thisRSID,:);
    TargetTT = GetTargetTT(ROOT,thisRSID,thisRegion,Params,0.03);

    thisTT_table = TT_table(TT_table.rat==thisRIDn & TT_table.session==thisSIDn,:);
    for t=1:size(thisTT_table,1)
        if ~ismember(thisTT_table.TT(t),TargetTT)
            thisTT_table.TT(t)=0;
        end
    end
    thisTT_table= thisTT_table(thisTT_table.TT~=0,:);
    [~,t] = max(thisTT_table.RippleBandMean);
    TargetTT_p = thisTT_table.TT(t);

    EEG = LoadEEGData(ROOT, thisRSID, TargetTT_p,Params,Params_Ripple);


    theseRip = RipplesTable(RipplesTable.rat==thisRIDn & RipplesTable.session==thisSIDn,:);

    if isempty(theseRip), continue; end

    for bid=1:size(BehavTable,1)-1

        try
        ITIs = BehavTable.xEnd(bid);
        ITIe = BehavTable.start(bid+1);

        thisRip = RipplesTable(RipplesTable.rat==thisRIDn & RipplesTable.session==thisSIDn&...
            RipplesTable.STtime>ITIs & RipplesTable.STtime<ITIe,:);

        Tori = theseRip.STtime(1) - theseRip.STindex(1) * unit;
        EPOCHst = floor((BehavTable.start(bid)-Tori) / unit );
        EPOCHed = ceil((BehavTable.reward(bid+1)-Tori )/ unit );
        thisEEG = EEG.(['TT' num2str(TargetTT_p)]).Raw(EPOCHst:EPOCHed);
        thisEEGr = EEG.(['TT' num2str(TargetTT_p)]).Filtered(EPOCHst:EPOCHed);

        EPOCHedT = (EPOCHed-EPOCHst)*unit;
        xlist = [0:5000:EPOCHed-EPOCHst]; xlistS={};
        for x=1:length(xlist),xlistS{x}=num2str(xlist(x)*unit); end
        %%
        clusters = UnitsTable(UnitsTable.rat==thisRIDn & UnitsTable.session==thisSIDn,:);
        cls_all = size(clusters,1);

        %%
        figure('position',[1927,82,1907,876]);
        subplot(5,1,1)
        plot(thisEEG,'color','k')
        line([1 1]*(BehavTable.xEnd(bid)-Tori) / unit - EPOCHst,[-2500 2500],'color','r')
        line([1 1]*(BehavTable.start(bid+1)-Tori) / unit - EPOCHst,[-2500 2500],'color','r')

        line([1 1]*(BehavTable.reward(bid)-Tori) / unit - EPOCHst,[-2500 2500],'color','b')
        line([1 1]*(BehavTable.reward(bid+1)-Tori) / unit - EPOCHst,[-2500 2500],'color','b')



        for rid=1:size(thisRip,1)
            x1 = (thisRip.STtime(rid)-Tori)/ unit - EPOCHst; x2 = (thisRip.EDtime(rid)-Tori)/ unit - EPOCHst;
            patch([x1 x2 x2 x1], [-1 -1 1 1]*2500,'y','facealpha',0.3,'edgealpha',0)
        end

        xticks(xlist)
        xticklabels(xlistS)
        xlim([0 xlist(end)])

        ylabel('Amplitude (uV)')
        %%
        subplot(5,1,2)
        plot(thisEEGr,'color','b')
        line([1 1]*(BehavTable.xEnd(bid)-Tori) / unit - EPOCHst,[-200 200],'color','r')
        line([1 1]*(BehavTable.start(bid+1)-Tori) / unit - EPOCHst,[-200 200],'color','r')

        line([1 1]*(BehavTable.reward(bid)-Tori) / unit - EPOCHst,[-200 200],'color','b')
        line([1 1]*(BehavTable.reward(bid+1)-Tori) / unit - EPOCHst,[-200 200],'color','b')

        for rid=1:size(thisRip,1)
            x1 = (thisRip.STtime(rid)-Tori)/ unit - EPOCHst; x2 = (thisRip.EDtime(rid)-Tori)/ unit - EPOCHst;
            patch([x1 x2 x2 x1], [-400 -400 400 400],'y','facealpha',0.3,'edgealpha',0)
        end

        xticks(xlist)
        xticklabels(xlistS)
         xlim([0 xlist(end)])
         ylabel('Amplitude (uV)')
%%
subplot(5,1,3)
thisPos = (Pos.t>=EPOCHst*unit+Tori & Pos.t<=EPOCHed*unit+Tori);
plot(max(Pos.y_linearized) - Pos.y_linearized(thisPos),'k','linewidth',2)
line([0 sum(thisPos)], [1 1]*max(Pos.y_linearized) -diverging_point,'color','b')
ylabel('Pos (cm)')
xlim([0 sum(thisPos)])

        %%
        subplot(5,1,[4 5]); hold on
        ylistS={};
        for un = 1:cls_all

            [thisRIDc,thisSIDc,thisTTIDc,thisCLIDc, thisFLIDc] = parsing_clusterID(clusters.ID{un},1);

            Spk = Spike.(['TT' (thisTTIDc)]).(['Unit' thisCLIDc]);

            thisSpks = Spk.t_spk(Spk.t_spk>=EPOCHst*unit+Tori & Spk.t_spk<=EPOCHed*unit+Tori);

            if clusters.(['Selectivity_' Sel])(un)>=1
                if (clusters.(['RDI_' Sel])(un))>0
                    c='r';
                else
                    c='b';
                end
            else
                c='k';
            end

            scatter(thisSpks-(EPOCHst*unit+Tori),ones(length(thisSpks),1)*un,c,'|')

            ylistS{un} = [(thisTTIDc) '-' thisCLIDc];
        end
        xlim([0 str2double(xlistS{end})])

        for rid=1:size(thisRip,1)
            x1 = thisRip.STtime(rid)-(EPOCHst*unit+Tori); x2 = thisRip.EDtime(rid)-(EPOCHst*unit+Tori);
            patch([x1 x2 x2 x1], [0 0 1 1]*(un+.5),'y','facealpha',0.3,'edgealpha',0)
        end

        yticks([1:un])
        yticklabels(ylistS)
        ylim([0 un+.5])

 ylabel('Unit')
        %%
        CxList = {'Zebra','Pebbles','Bamboo','Mountains'};
        CrList={'C','W'};
        cx1 = CxList{BehavTable.context(bid)}; cx2 = CxList{BehavTable.context(bid+1)}; 
        cr1 = CrList{BehavTable.correctness(bid)}; cr2 = CrList{BehavTable.correctness(bid+1)}; 
        sgtitle(['Trial ' BehavTable.TrialID{bid} ' (' cx1 '(' cr1 ') -' cx2 '(' cr2 '))'])


        ROOT.Fig_en = [ROOT.Fig '\All\' thisRegion '\' Sel];
        if ~exist(ROOT.Fig_en), mkdir(ROOT.Fig_en); end
        saveas(gca,[ROOT.Fig_en '\'  BehavTable.TrialID{bid} '.svg'])
        saveas(gca,[ROOT.Fig_en '\'  BehavTable.TrialID{bid} '.png'])

        if ~isempty(thisRip)
                ROOT.Fig_en = [ROOT.Fig '\Rip\' thisRegion '\' Sel];
        if ~exist(ROOT.Fig_en), mkdir(ROOT.Fig_en); end
        saveas(gca,[ROOT.Fig_en '\'  BehavTable.TrialID{bid} '.svg'])
        saveas(gca,[ROOT.Fig_en '\'  BehavTable.TrialID{bid} '.png'])
        end
        
        close all
        catch
            close all
        end
    end

end
