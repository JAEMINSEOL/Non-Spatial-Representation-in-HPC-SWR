Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data (Blur)'];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R5'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R4'];
ROOT.Units = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];
if ~exist(ROOT.Fig3), mkdir(ROOT.Fig3); end
if ~exist(ROOT.Rip), mkdir(ROOT.Rip); end

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
randN=5000;
%%
RipplesTable_p = RipplesTable(RipplesTable.nPCs>4,:);
for rid=1:size(RipplesTable_p,1)
    RipID = RipplesTable_p.ID{rid};
    thisReact = ReactTable(strcmp(ReactTable.RippleID,RipID),:);
    thisUnits=table;
    for u=1:size(thisReact,1)
        thisUnits =[thisUnits;UnitsTable_A(find(cellfun(Params.cellfind(thisReact.UnitID(u)),UnitsTable_A.ID)),:)];
    end

    thisUnits_All = thisUnits(~isnan(thisUnits.RDI_All),:);
    thisUnits_No = thisUnits(~isnan(thisUnits.RDI_No),:);
    thisUnits_Lo = thisUnits(~isnan(thisUnits.RDI_Lo),:);
    thisUnits_Hi = thisUnits(~isnan(thisUnits.RDI_Hi),:);

    thisPool = UnitsTable(UnitsTable.rat==RipplesTable_p.rat(rid) & UnitsTable.session==RipplesTable_p.session(rid),:);

    thisPool_All =  thisPool(~isnan(thisPool.RDI_All),:);
    thisPool_No =  thisPool(~isnan(thisPool.RDI_No),:);
    thisPool_Lo =  thisPool(~isnan(thisPool.RDI_Lo),:);
    thisPool_Hi =  thisPool(~isnan(thisPool.RDI_Hi),:);

    RDI_All.dist = table;
    RDI_No.dist=table;
    RDI_Lo.dist=table;
    RDI_Hi.dist=table;

    RDI_All.dist.mean=nan;
    RDI_No.dist.mean=nan;
    RDI_Lo.dist.mean=nan;
    RDI_Hi.dist.mean=nan;

    RDI_All.dist.median=nan;
    RDI_No.dist.median=nan;
    RDI_Lo.dist.median=nan;
    RDI_Hi.dist.median=nan;

    try
        for i=1:randN
            samp = datasample(thisPool_All,size(thisUnits_All,1),1);
            if ~isempty(samp), RDI_All.dist.mean(i) = nanmean(samp.RDI_All); end
            samp = datasample(thisPool_No,size(thisUnits_No,1),1);
            if ~isempty(samp), RDI_No.dist.mean(i) = nanmean(samp.RDI_No); end
            samp = datasample(thisPool_Lo,size(thisUnits_Lo,1),1);
            if ~isempty(samp), RDI_Lo.dist.mean(i) = nanmean(samp.RDI_Lo); end
            samp = datasample(thisPool_Hi,size(thisUnits_Hi,1),1);
            if ~isempty(samp), RDI_Hi.dist.mean(i) = nanmean(samp.RDI_Hi); end

            if ~isempty(samp), RDI_All.dist.median(i) = nanmedian(samp.RDI_All); end
            if ~isempty(samp), RDI_No.dist.median(i) = nanmedian(samp.RDI_No); end
            if ~isempty(samp), RDI_Lo.dist.median(i) = nanmedian(samp.RDI_Lo); end
            if ~isempty(samp), RDI_Hi.dist.median(i) = nanmedian(samp.RDI_Hi); end
        end
    catch
        disp([RipID ' perm fail'])
    end
    RDI_All.act_mean = nanmean(thisUnits.RDI_All);
    RDI_No.act_mean = nanmean(thisUnits.RDI_No);
    RDI_Lo.act_mean = nanmean(thisUnits.RDI_Lo);
    RDI_Hi.act_mean = nanmean(thisUnits.RDI_Hi);

    RDI_All.act_median = nanmedian(thisUnits.RDI_All);
    RDI_No.act_median = nanmedian(thisUnits.RDI_No);
    RDI_Lo.act_median = nanmedian(thisUnits.RDI_Lo);
    RDI_Hi.act_median = nanmedian(thisUnits.RDI_Hi);

    RDI_All.p_mean = min(sum(RDI_All.act_mean>RDI_All.dist.mean),sum(RDI_All.act_mean<RDI_All.dist.mean))/sum(~isnan(RDI_All.dist.mean));
    RDI_No.p_mean = min(sum(RDI_No.act_mean>RDI_No.dist.mean),sum(RDI_No.act_mean<RDI_No.dist.mean))/sum(~isnan(RDI_No.dist.mean));
    RDI_Lo.p_mean = min(sum(RDI_Lo.act_mean>RDI_Lo.dist.mean),sum(RDI_Lo.act_mean<RDI_Lo.dist.mean))/sum(~isnan(RDI_Lo.dist.mean));
    RDI_Hi.p_mean = min(sum(RDI_Hi.act_mean>RDI_Hi.dist.mean),sum(RDI_Hi.act_mean<RDI_Hi.dist.mean))/sum(~isnan(RDI_Hi.dist.mean));

    RDI_All.p_median = min(sum(RDI_All.act_median>RDI_All.dist.median),sum(RDI_All.act_median<RDI_All.dist.median))/...
        sum(~isnan(RDI_All.dist.median));
    RDI_No.p_median = min(sum(RDI_No.act_median>RDI_No.dist.median),sum(RDI_No.act_median<RDI_No.dist.median))/...
        sum(~isnan(RDI_No.dist.median));
    RDI_Lo.p_median = min(sum(RDI_Lo.act_median>RDI_Lo.dist.median),sum(RDI_Lo.act_median<RDI_Lo.dist.median))/...
        sum(~isnan(RDI_Lo.dist.median));
    RDI_Hi.p_median = min(sum(RDI_Hi.act_median>RDI_Hi.dist.median),sum(RDI_Hi.act_median<RDI_Hi.dist.median))/...
        sum(~isnan(RDI_Hi.dist.median));

    save([ROOT.Rip ['\R-' RipID '.mat']],'thisUnits','RDI_All','RDI_No','RDI_Lo','RDI_Hi')
    disp([RipID ' is finished!'])

    RipplesTable_p.pRDI_All(rid) = RDI_All.p_mean;
    RipplesTable_p.pRDI_No(rid) = RDI_No.p_mean;
    RipplesTable_p.pRDI_Lo(rid) = RDI_Lo.p_mean;
    RipplesTable_p.pRDI_Hi(rid) = RDI_Hi.p_mean;

end
 suff = '_RDIs';
writetable(RipplesTable_p,[ROOT.Save '\RipplesTable_' thisRegion suff '.xlsx'],'writemode','overwrite')