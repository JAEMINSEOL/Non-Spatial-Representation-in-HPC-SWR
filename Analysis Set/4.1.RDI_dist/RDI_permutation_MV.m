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

parpool('local',8); % Change 4 to the number of workers you want
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'SUB';
thisRegion2 = 'SUB_field';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion2 '.xlsx']);

UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
randN=5000;
%%
col_name_A = string(UnitsTable_A.Properties.VariableNames(24:32));
UnitsTable_B(1,28:36) = array2table(nan(1,9));
UnitsTable_B = renamevars(UnitsTable_B, UnitsTable_B.Properties.VariableNames(28:36),col_name_A);
for uid=1:size(UnitsTable_B,1)
    tar = find(strncmp(UnitsTable_A.ID,UnitsTable_B.ID(uid),12));
    UnitsTable_B(uid,28:36) = UnitsTable_A(tar,24:32);
end
    
UnitsTable = UnitsTable_B;
%%
% RipplesTable_p = RipplesTable(RipplesTable.nPCs>2,:);
RipplesTable_p = RipplesTable;
for r=1:size(RipplesTable_p,1)
    RipID = RipplesTable_p.ID{r};
    thisReact = ReactTable(strcmp(ReactTable.RippleID,RipID),:);
    thisUnits=table;
    for u=1:size(thisReact,1)
        thisUnits =[thisUnits;UnitsTable(find(cellfun(Params.cellfind(thisReact.UnitID(u)),UnitsTable.ID)),:)];
    end


    if size(thisUnits,1)==0

        RipplesTable_p.pRDI_L(r)=nan;
        RipplesTable_p.pRDI_R(r)=nan;
        RipplesTable_p.pRDI_C(r)=nan;

        RipplesTable_p.pRDI_L_scuv(r)=nan;
        RipplesTable_p.pRDI_R_scuv(r)=nan;
        RipplesTable_p.pRDI_C_scuv(r)=nan;
        continue;
    end



    thisUnits_UV_L = thisUnits((~isnan(thisUnits.RDI_LScene) & ~thisUnits.MultiVar_L),:);
    thisUnits_UV_R = thisUnits((~isnan(thisUnits.RDI_RScene) & ~thisUnits.MultiVar_R),:);
    thisUnits_UV_C = thisUnits((~isnan(thisUnits.RDI_LR) & ~thisUnits.MultiVar_C),:);
    thisUnits_UV_SC = thisUnits((~isnan(thisUnits.RDI_LR) & ~thisUnits.MultiVar_SC),:);

    thisPool = UnitsTable(UnitsTable.rat==RipplesTable_p.rat(r) & UnitsTable.session==RipplesTable_p.session(r),:);

    thisPool_UV_L =  thisPool((~isnan(thisPool.RDI_LScene)& ~thisPool.MultiVar_L),:);
    thisPool_UV_R = thisPool((~isnan(thisPool.RDI_RScene)& ~thisPool.MultiVar_R),:);
    thisPool_UV_C =  thisPool((~isnan(thisPool.RDI_LR)& ~thisPool.MultiVar_C),:);
    thisPool_UV_SC =  thisPool((~isnan(thisPool.RDI_LR)& ~thisPool.MultiVar_SC),:);


    RDI_L.dist=table;
    RDI_R.dist=table;
    RDI_C.dist=table;

    temp1=nan; temp2=nan; temp3=nan; temp4=nan; temp5=nan; temp6=nan;
    temp1_scuv=nan; temp2_scuv=nan; temp3_scuv=nan; temp4_scuv=nan; temp5_scuv=nan; temp6_scuv=nan;

    tic
    try
        parfor i=1:randN
            samp = datasample(thisPool_UV_L,size(thisUnits_UV_L,1),1);
            if ~isempty(samp), temp1(i) = nanmean(samp.RDI_LScene); temp2(i) = nanmedian(samp.RDI_LScene); end
            samp = datasample(thisPool_UV_R,size(thisUnits_UV_R,1),1);
            if ~isempty(samp), temp3(i) = nanmean(samp.RDI_RScene); temp4(i) = nanmedian(samp.RDI_RScene); end
            samp = datasample(thisPool_UV_C,size(thisUnits_UV_C,1),1);
            if ~isempty(samp), temp5(i) = nanmean(samp.RDI_LR); temp6(i) = nanmedian(samp.RDI_LR); end

            samp = datasample(thisPool_UV_SC,size(thisUnits_UV_L,1),1);
            if ~isempty(samp),  temp1_scuv(i) = nanmean(samp.RDI_LScene);  temp2_scuv(i) = nanmedian(samp.RDI_LScene);  end
            if ~isempty(samp),  temp3_scuv(i) = nanmean(samp.RDI_RScene); temp4_scuv(i) = nanmedian(samp.RDI_RScene); end
            if ~isempty(samp), temp5_scuv(i) = nanmean(samp.RDI_LR); temp6_scuv(i) = nanmedian(samp.RDI_LR); end
        end
    catch
        disp([RipID ' perm fail'])
    end
    toc

    RDI_L.dist.mean = temp1;
    RDI_L.dist.median = temp2;
    RDI_R.dist.mean = temp3;
    RDI_R.dist.median = temp4;
    RDI_C.dist.mean = temp5;
    RDI_C.dist.median = temp6;

    RDI_L.dist.mean_scuv = temp1_scuv;
    RDI_L.dist.median_scuv = temp2_scuv;
    RDI_R.dist.mean_scuv = temp3_scuv;
    RDI_R.dist.median_scuv = temp4_scuv;
    RDI_C.dist.mean_scuv = temp5_scuv;
    RDI_C.dist.median_scuv = temp6_scuv;


    RDI_L.act_mean = nanmean(thisUnits_UV_L.RDI_LScene);
    RDI_R.act_mean = nanmean(thisUnits_UV_R.RDI_RScene);
    RDI_C.act_mean = nanmean(thisUnits_UV_C.RDI_LR);
    RDI_L.act_mean_scuv = nanmean(thisUnits_UV_SC.RDI_LScene);
    RDI_R.act_mean_scuv = nanmean(thisUnits_UV_SC.RDI_RScene);
    RDI_C.act_mean_scuv = nanmean(thisUnits_UV_SC.RDI_LR);

    RDI_L.act_median = nanmedian(thisUnits_UV_L.RDI_LScene);
    RDI_R.act_median = nanmedian(thisUnits_UV_R.RDI_RScene);
    RDI_C.act_median = nanmedian(thisUnits_UV_C.RDI_LR);
    RDI_L.act_median_scuv = nanmedian(thisUnits_UV_SC.RDI_LScene);
    RDI_R.act_median_scuv = nanmedian(thisUnits_UV_SC.RDI_RScene);
    RDI_C.act_median_scuv = nanmedian(thisUnits_UV_SC.RDI_LR);

    RDI_L.p_mean = min(sum(RDI_L.act_mean>RDI_L.dist.mean),sum(RDI_L.act_mean<RDI_L.dist.mean))/sum(~isnan(RDI_L.dist.mean));
    RDI_R.p_mean = min(sum(RDI_R.act_mean>RDI_R.dist.mean),sum(RDI_R.act_mean<RDI_R.dist.mean))/sum(~isnan(RDI_R.dist.mean));
    RDI_C.p_mean = min(sum(RDI_C.act_mean>RDI_C.dist.mean),sum(RDI_C.act_mean<RDI_C.dist.mean))/sum(~isnan(RDI_C.dist.mean));

    RDI_L.p_mean_scuv = min(sum(RDI_L.act_mean_scuv>RDI_L.dist.mean_scuv),sum(RDI_L.act_mean_scuv<RDI_L.dist.mean_scuv))/sum(~isnan(RDI_L.dist.mean_scuv));
    RDI_R.p_mean_scuv = min(sum(RDI_R.act_mean_scuv>RDI_R.dist.mean_scuv),sum(RDI_R.act_mean_scuv<RDI_R.dist.mean_scuv))/sum(~isnan(RDI_R.dist.mean_scuv));
    RDI_C.p_mean_scuv = min(sum(RDI_C.act_mean_scuv>RDI_C.dist.mean_scuv),sum(RDI_C.act_mean_scuv<RDI_C.dist.mean_scuv))/sum(~isnan(RDI_C.dist.mean_scuv));

    RDI_L.p_median = min(sum(RDI_L.act_median>RDI_L.dist.median),sum(RDI_L.act_median<RDI_L.dist.median))/...
        sum(~isnan(RDI_L.dist.median));
    RDI_R.p_median = min(sum(RDI_R.act_median>RDI_R.dist.median),sum(RDI_R.act_median<RDI_R.dist.median))/...
        sum(~isnan(RDI_R.dist.median));
    RDI_C.p_median = min(sum(RDI_C.act_median>RDI_C.dist.median),sum(RDI_C.act_median<RDI_C.dist.median))/...
        sum(~isnan(RDI_C.dist.median));

    RDI_L.p_median_scuv = min(sum(RDI_L.act_median_scuv>RDI_L.dist.median_scuv),sum(RDI_L.act_median_scuv<RDI_L.dist.median_scuv))/...
        sum(~isnan(RDI_L.dist.median_scuv));
    RDI_R.p_median_scuv = min(sum(RDI_R.act_median_scuv>RDI_R.dist.median_scuv),sum(RDI_R.act_median_scuv<RDI_R.dist.median_scuv))/...
        sum(~isnan(RDI_R.dist.median_scuv));
    RDI_C.p_median_scuv = min(sum(RDI_C.act_median_scuv>RDI_C.dist.median_scuv),sum(RDI_C.act_median_scuv<RDI_C.dist.median_scuv))/...
        sum(~isnan(RDI_C.dist.median_scuv));

    save([ROOT.Rip ['\' thisRegion2 '-' RipID '_UV.mat']],'thisUnits','RDI_L','RDI_R','RDI_C')
    disp([RipID ' is finished!'])


    RipplesTable_p.pRDI_L_UV(r) = RDI_L.p_mean;
    RipplesTable_p.pRDI_R_UV(r) = RDI_R.p_mean;
    RipplesTable_p.pRDI_C_UV(r) = RDI_C.p_mean;

        RipplesTable_p.pRDI_L_UV_scuv(r) = RDI_L.p_mean_scuv;
    RipplesTable_p.pRDI_R_UV_scuv(r) = RDI_R.p_mean_scuv;
    RipplesTable_p.pRDI_C_UV_scuv(r) = RDI_C.p_mean_scuv;


    if size(thisUnits_UV_L,1)<5, RipplesTable_p.pRDI_L_UV(r)=nan; end
    if size(thisUnits_UV_R,1)<5, RipplesTable_p.pRDI_R_UV(r)=nan; end
    if size(thisUnits_UV_C,1)<5, RipplesTable_p.pRDI_C_UV(r)=nan; end

    if size(thisUnits_UV_SC,1)<5, RipplesTable_p.pRDI_L_UV_scuv(r)=nan; RipplesTable_p.pRDI_R_UV_scuv(r)=nan; RipplesTable_p.pRDI_C_UV_scuv(r)=nan; end

end
suff = '_RDIs_UV';
writetable(RipplesTable_p,[ROOT.Save '\RipplesTable_' thisRegion2 suff '.xlsx'],'writemode','replacefile')
delete(gcp('nocreate'));