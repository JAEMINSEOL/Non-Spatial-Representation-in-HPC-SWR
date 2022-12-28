Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R5'];
ROOT.Fig3 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R4'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
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
RipplesTable_p = RipplesTable(RipplesTable.nPCs>2,:);
for r=1:size(RipplesTable_p,1)
    RipID = RipplesTable_p.ID{r};
    thisReact = ReactTable(strcmp(ReactTable.RippleID,RipID),:);
    thisUnits=table;
    for u=1:size(thisReact,1)
    thisUnits =[thisUnits;UnitsTable_A(find(cellfun(Params.cellfind(thisReact.UnitID(u)),UnitsTable_A.ID)),:)];
    end

    thisUnits_L = thisUnits(~isnan(thisUnits.RDI_LScene),:);
    thisUnits_R = thisUnits(~isnan(thisUnits.RDI_RScene),:);
    thisUnits_C = thisUnits(~isnan(thisUnits.RDI_LR),:);

    thisPool = UnitsTable(UnitsTable.rat==RipplesTable_p.rat(r) & UnitsTable.session==RipplesTable_p.session(r),:);

         thisPool_L =  thisPool(~isnan(thisPool.RDI_LScene),:);
    thisPool_R =  thisPool(~isnan(thisPool.RDI_RScene),:);
    thisPool_C =  thisPool(~isnan(thisPool.RDI_LR),:);


    RDI_L.dist=table;
    RDI_R.dist=table;
    RDI_C.dist=table;

    try
for i=1:randN
    samp = datasample(thisPool_L,size(thisUnits_L,1),1);
    RDI_L.dist.mean(i) = nanmean(samp.RDI_LScene);
    samp = datasample(thisPool_R,size(thisUnits_R,1),1);
    RDI_R.dist.mean(i) = nanmean(samp.RDI_RScene);
    samp = datasample(thisPool_C,size(thisUnits_C,1),1);
    RDI_C.dist.mean(i) = nanmean(samp.RDI_LR);

        RDI_L.dist.median(i) = nanmedian(samp.RDI_LScene);
    RDI_R.dist.median(i) = nanmedian(samp.RDI_RScene);
    RDI_C.dist.median(i) = nanmedian(samp.RDI_LR);
end
    catch
        disp([RipID ' perm fail'])
    end
    
    RDI_L.act_mean = nanmean(thisUnits.RDI_LScene);
    RDI_R.act_mean = nanmean(thisUnits.RDI_RScene);
    RDI_C.act_mean = nanmean(thisUnits.RDI_LR);

            RDI_L.act_median = nanmedian(thisUnits.RDI_LScene);
    RDI_R.act_median = nanmedian(thisUnits.RDI_RScene);
    RDI_C.act_median = nanmedian(thisUnits.RDI_LR);

    RDI_L.p_mean = min(sum(RDI_L.act_mean>RDI_L.dist.mean),sum(RDI_L.act_mean<RDI_L.dist.mean))/sum(~isnan(RDI_L.dist.mean));
    RDI_R.p_mean = min(sum(RDI_R.act_mean>RDI_R.dist.mean),sum(RDI_R.act_mean<RDI_R.dist.mean))/sum(~isnan(RDI_R.dist.mean));
    RDI_C.p_mean = min(sum(RDI_C.act_mean>RDI_C.dist.mean),sum(RDI_C.act_mean<RDI_C.dist.mean))/sum(~isnan(RDI_C.dist.mean));

        RDI_L.p_median = min(sum(RDI_L.act_median>RDI_L.dist.median),sum(RDI_L.act_median<RDI_L.dist.median))/...
            sum(~isnan(RDI_L.dist.median));
    RDI_R.p_median = min(sum(RDI_R.act_median>RDI_R.dist.median),sum(RDI_R.act_median<RDI_R.dist.median))/...
        sum(~isnan(RDI_R.dist.median));
    RDI_C.p_median = min(sum(RDI_C.act_median>RDI_C.dist.median),sum(RDI_C.act_median<RDI_C.dist.median))/...
        sum(~isnan(RDI_C.dist.median));

    save([ROOT.Rip ['\R-' RipID '.mat']],'thisUnits','RDI_L','RDI_R','RDI_C')
    disp([RipID ' is finished!'])
end

    
