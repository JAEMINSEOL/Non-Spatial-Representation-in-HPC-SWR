Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
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
thisRegion2 = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '.xlsx']);
% UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion2 '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']); 
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
randN=5000;
%
% RipplesTable_p = RipplesTable(RipplesTable.nPCs>2,:);
RipplesTable_p = RipplesTable;
for r=1:size(RipplesTable_p,1)
    RipID = RipplesTable_p.ID{r};
    thisReact = ReactTable(strcmp(ReactTable.RippleID,RipID),:);
    thisUnits=table;
    for u=1:size(thisReact,1)
        thisUnits =[thisUnits;UnitsTable_A(find(cellfun(Params.cellfind(thisReact.UnitID(u)),UnitsTable_A.ID)),:)];
    end

%     thisUnits_SF = thisUnits(thisUnits.NumField<2,:);
thisUnits_SF = thisUnits;
    thisUnits = thisUnits_SF;

    if size(thisUnits_SF,1)==0

        RipplesTable_p.pRDI_L(r)=nan;
        RipplesTable_p.pRDI_R(r)=nan;
        RipplesTable_p.pRDI_C(r)=nan;
        continue; 
    end



    thisUnits_L = thisUnits_SF(~isnan(thisUnits_SF.RDI_LScene) & abs(thisUnits_SF.RDI_LScene)>0.1,:);
    thisUnits_R = thisUnits_SF(~isnan(thisUnits_SF.RDI_RScene) & abs(thisUnits_SF.RDI_RScene)>0.1,:);
    thisUnits_C = thisUnits_SF(~isnan(thisUnits_SF.RDI_LR) & abs(thisUnits_SF.RDI_LR)>0.1,:);

    thisPool = UnitsTable(UnitsTable.rat==RipplesTable_p.rat(r) & UnitsTable.session==RipplesTable_p.session(r),:);

    thisPool_SF = thisPool;
%     thisPool_SF = thisPool(thisPool.NumField<2,:);

    thisPool_L =  thisPool_SF(~isnan(thisPool_SF.RDI_LScene) & abs(thisPool_SF.RDI_LScene)>0.1,:);
    thisPool_R =  thisPool_SF(~isnan(thisPool_SF.RDI_RScene) & abs(thisPool_SF.RDI_RScene)>0.1,:);
    thisPool_C =  thisPool_SF(~isnan(thisPool_SF.RDI_LR) & abs(thisPool_SF.RDI_LR)>0.1,:);

    thisPool_L =  thisPool_SF(~isnan(thisPool_SF.RDI_LScene),:);
    thisPool_R =  thisPool_SF(~isnan(thisPool_SF.RDI_RScene),:);
    thisPool_C =  thisPool_SF(~isnan(thisPool_SF.RDI_LR),:);


    RDI_L.dist=table;
    RDI_R.dist=table;
    RDI_C.dist=table;
      RDI_L.dist.mean=nan;
    RDI_R.dist.mean=nan;
    RDI_C.dist.mean=nan;
    RDI_L.dist.median=nan;
    RDI_R.dist.median=nan;
    RDI_C.dist.median=nan;



    temp1=nan; temp2=nan; temp3=nan; temp4=nan; temp5=nan; temp6=nan;
        try
        parfor i=1:randN
            samp = datasample(thisPool_L,size(thisUnits_L,1),1);
            if ~isempty(samp), temp1(i) = nanmean(samp.RDI_LScene); temp2(i) = nanmedian(samp.RDI_LScene); end
            samp = datasample(thisPool_R,size(thisUnits_R,1),1);
            if ~isempty(samp), temp3(i) = nanmean(samp.RDI_RScene); temp4(i) = nanmedian(samp.RDI_RScene); end
            samp = datasample(thisPool_C,size(thisUnits_C,1),1);
            if ~isempty(samp), temp5(i) = nanmean(samp.RDI_LR); temp6(i) = nanmedian(samp.RDI_LR); end

        end
    catch
        disp([RipID ' perm fail'])
        end

            RDI_L.dist.mean = temp1;
    RDI_L.dist.median = temp2;
    RDI_R.dist.mean = temp3;
    RDI_R.dist.median = temp4;
    RDI_C.dist.mean = temp5;
    RDI_C.dist.median = temp6;


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

    save([ROOT.Rip ['\' thisRegion2 '-' RipID '_SF.mat']],'thisUnits','RDI_L','RDI_R','RDI_C')
    disp([RipID ' is finished!'])

   RipplesTable_p.pRDI_L(r) = RDI_L.p_mean;
    RipplesTable_p.pRDI_R(r) = RDI_R.p_mean;
    RipplesTable_p.pRDI_C(r) = RDI_C.p_mean;

        if size(thisUnits_L,1)<5, RipplesTable_p.pRDI_L(r)=nan; end
    if size(thisUnits_R,1)<5, RipplesTable_p.pRDI_R(r)=nan; end
    if size(thisUnits_C,1)<5, RipplesTable_p.pRDI_C(r)=nan; end
%     RipplesTable_p.nRDI_hetero(r) = nanmean(thisUnits.MultiVar);
%     RipplesTable_p.mRDI_hetero(r) = nanmean(thisUnits.RDI_hetero);

end
 suff = '_RDIs_Sig';
writetable(RipplesTable_p,[ROOT.Save '\RipplesTable_' thisRegion2 suff '.xlsx'],'writemode','replacefile')

RipplesTable_p = readtable([ROOT.Save '\RipplesTable_' thisRegion2 suff '.xlsx']);
ns = isnan(RipplesTable_p.pRDI_L) & isnan(RipplesTable_p.pRDI_R) & isnan(RipplesTable_p.pRDI_C);
RipplesTable_p = RipplesTable_p(~ns,:);
    nsid = nanmin([RipplesTable_p.pRDI_L,RipplesTable_p.pRDI_R,RipplesTable_p.pRDI_C],[],2)<0.05;
    sum(RipplesTable_p.DecodingP_all<0.05 & nsid)

    RipplesTable_p = RipplesTable_p(RipplesTable_p.correctness,:);
