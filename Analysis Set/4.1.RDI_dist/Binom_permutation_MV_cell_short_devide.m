addpath('D:\HPC-SWR project\Analysis Program')
Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R4'];
ROOT.Units = [ROOT.Save '\units_mat\U2'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];
if ~exist(ROOT.Fig3), mkdir(ROOT.Fig3); end
if ~exist(ROOT.Rip), mkdir(ROOT.Rip); end

% parpool('local',8); % Change 4 to the number of workers you want
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion0 = 'CA1';
thisRegion = 'CA1';
thisRegion2 = 'CA1_field';

% thisRegion0 = 'SUB';
% thisRegion = 'SUB';
% thisRegion2 = 'SUB_field';


% RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis' '.xlsx']);
ReactTableA = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);
ReactTableB = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion2 '.xlsx']);


UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);

TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
randN=2000;

UnitsTable = UnitsTable_B;

%%
% RipplesTable_p = RipplesTable(RipplesTable.nPCs>2,:);
RipplesTable_p = RipplesTable;
for rid=1:size(RipplesTable_p,1)
    RipID = RipplesTable_p.ID{rid};
    thisReact_A = ReactTableA(strcmp(ReactTableA.RippleID,RipID),:);
    thisReact_A = unique(thisReact_A(:,1:2));

        thisReact_B = ReactTableB(strcmp(ReactTableB.RippleID,RipID),:);
    thisReact_B = unique(thisReact_B(:,1:2));

    for uid=1:size(thisReact_A,1)
        try
        idx=find(strcmp(thisReact_A.UnitID{uid},UnitsTable_A.ID));
        thisReact_A.sel_L(uid)=UnitsTable_A.Selectivity_LScene(idx);
         thisReact_A.sel_R(uid)=UnitsTable_A.Selectivity_RScene(idx);
 thisReact_A.sel_C(uid)=UnitsTable_A.Selectivity_LR(idx);
        catch
        end
    end

    for uid=1:size(thisReact_B,1)
        try
        idx=find(strcmp(thisReact_B.UnitID{uid},UnitsTable_B.ID));
        thisReact_B.sel_L(uid)=UnitsTable_B.RDI_LScene(idx);
         thisReact_B.sel_R(uid)=UnitsTable_B.RDI_RScene(idx);
 thisReact_B.sel_C(uid)=UnitsTable_B.RDI_LR(idx);
        catch
        end
    end


try
%     if exist([ROOT.Rip ['\' thisRegion0 '-' RipID '_UV.mat']])

%         load([ROOT.Rip ['\' thisRegion0 '-' RipID '_UV.mat']],'thisUnits','RDI_L','RDI_R','RDI_C')
        %     disp([RipID ' is finished!'])

    thisUnits=table;
    for u=1:size(thisReact_B,1)
        thisUnits =[thisUnits;UnitsTable(find(cellfun(Params.cellfind(thisReact_B.UnitID(u)),UnitsTable.ID)),:)];
    end

%         RipplesTable_p.mRDI_L_UV(rid) = RDI_L.act_mean;
%         RipplesTable_p.mRDI_R_UV(rid) = RDI_R.act_mean;
%         RipplesTable_p.mRDI_C_UV(rid) = RDI_C.act_mean;

m = median(UnitsTable_B.RDI_LScene(abs(UnitsTable_B.RDI_LScene)>=0.1));
p1 = sum(UnitsTable_B.RDI_LScene>=0.1)/sum(abs(UnitsTable_B.RDI_LScene)>=0.1); p2 = sum(UnitsTable_B.RDI_LScene<=-0.1)/sum(abs(UnitsTable_B.RDI_LScene)>=0.1);
if sum(abs(thisUnits.RDI_LScene)>=0.1)<1
    RipplesTable_p.WilRDI_L_UV(rid) = nan; 
    RipplesTable_p.MRDI_L_UV(rid) = nan;
    RipplesTable_p.pBinom_L_UV(rid) = nan;
else
[p,h,w] = signrank(thisUnits.RDI_LScene(abs(thisUnits.RDI_LScene)>=0.1),m);
RipplesTable_p.WilRDI_L_UV(rid) = p;
RipplesTable_p.MRDI_L_UV(rid) = median(thisUnits.RDI_LScene(abs(thisUnits.RDI_LScene)>=0.1));
b1=myBinomTest(sum(thisUnits.RDI_LScene>=0.1),sum(abs(thisUnits.RDI_LScene)>=0.1),p1,'one');
b2= myBinomTest(sum(thisUnits.RDI_LScene<=-0.1),sum(abs(thisUnits.RDI_LScene)>=0.1),p2,'one');
RipplesTable_p.pBinom_L_UV(rid) = min([b1,b2]);
end

m = median(UnitsTable_B.RDI_RScene(abs(UnitsTable_B.RDI_RScene)>=0.1));
p1 = sum(UnitsTable_B.RDI_RScene>=0.1)/sum(abs(UnitsTable_B.RDI_RScene)>=0.1); p2 = sum(UnitsTable_B.RDI_RScene<=-0.1)/sum(abs(UnitsTable_B.RDI_RScene)>=0.1);
if sum(abs(thisUnits.RDI_RScene)>=0.1)<1
    RipplesTable_p.WilRDI_R_UV(rid) = nan;
    RipplesTable_p.MRDI_R_UV(rid) = nan;
     RipplesTable_p.pBinom_R_UV(rid) = nan;
else
[p,h,w] = signrank(thisUnits.RDI_RScene(abs(thisUnits.RDI_RScene)>=0.1)-m);
RipplesTable_p.WilRDI_R_UV(rid) = p;
RipplesTable_p.MRDI_R_UV(rid) = median(thisUnits.RDI_RScene(abs(thisUnits.RDI_RScene)>=0.1));
b1=myBinomTest(sum(thisUnits.RDI_RScene>=0.1),sum(abs(thisUnits.RDI_RScene)>=0.1),p1,'one');
b2= myBinomTest(sum(thisUnits.RDI_RScene<=-0.1),sum(abs(thisUnits.RDI_RScene)>=0.1),p2,'one');
RipplesTable_p.pBinom_R_UV(rid) = min([b1,b2]);
end

m = median(UnitsTable_B.RDI_LR(abs(UnitsTable_B.RDI_LR)>=0.1));
p1 = sum(UnitsTable_B.RDI_LR>=0.1)/sum(abs(UnitsTable_B.RDI_LR)>=0.1); p2 = sum(UnitsTable_B.RDI_LR<=-0.1)/sum(abs(UnitsTable_B.RDI_LR)>=0.1);
if sum(abs(thisUnits.RDI_LR)>=0.1)<3
    RipplesTable_p.WilRDI_C_UV(rid) = nan;
    RipplesTable_p.MRDI_C_UV(rid) = nan;
     RipplesTable_p.pBinom_C_UV(rid) = nan;
else
[p,h,w] = signrank(thisUnits.RDI_LR(abs(thisUnits.RDI_LR)>=0.1)-m);
RipplesTable_p.WilRDI_C_UV(rid) = p;
RipplesTable_p.MRDI_C_UV(rid) = median(thisUnits.RDI_LR(abs(thisUnits.RDI_LR)>=0.1));
b1=myBinomTest(sum(thisUnits.RDI_LR>=0.1),sum(abs(thisUnits.RDI_LR)>=0.1),p1,'one');
b2= myBinomTest(sum(thisUnits.RDI_LR<=-0.1),sum(abs(thisUnits.RDI_LR)>=0.1),p2,'one');
RipplesTable_p.pBinom_C_UV(rid) = min([b1,b2]);
end

%  sum(nanmin([RipplesTable_p.WilRDI_L_UV, RipplesTable_p.WilRDI_R_UV RipplesTable_p.WilRDI_C_UV],[],2)<0.05)
%  sum(nanmin([RipplesTable_p.pBinom_L_UV, RipplesTable_p.pBinom_R_UV RipplesTable_p.pBinom_C_UV],[],2)<0.05)
%         RipplesTable_p.mRDI_L_UV(rid) = mean(thisUnits.RDI_LScene(abs(thisUnits.RDI_LScene)>=0.1));
%         RipplesTable_p.mRDI_R_UV(rid) = mean(thisUnits.RDI_RScene(abs(thisUnits.RDI_RScene)>=0.1));
%         RipplesTable_p.mRDI_C_UV(rid) = mean(thisUnits.RDI_LR(abs(thisUnits.RDI_LR)>=0.1));
% 
% 
%        RipplesTable_p.snr_RDI_L_UV(rid) = cal_snr(thisUnits.RDI_LScene);
%        RipplesTable_p.snr_RDI_R_UV(rid) = cal_snr(thisUnits.RDI_RScene); 
%        RipplesTable_p.snr_RDI_C_UV(rid) = cal_snr(thisUnits.RDI_LR); 
% 
% 
% 
%         if isnan(RipplesTable_p.mRDI_L_UV(rid)), RipplesTable_p.mRDI_L_UV(rid) = nanmean(thisUnits.RDI_LScene); end
%         if isnan(RipplesTable_p.mRDI_R_UV(rid)), RipplesTable_p.mRDI_R_UV(rid) = nanmean(thisUnits.RDI_RScene); end
%         if isnan(RipplesTable_p.mRDI_C_UV(rid)), RipplesTable_p.mRDI_C_UV(rid) = nanmean(thisUnits.RDI_LR); end
% 
%         RipplesTable_p.n_sf_L(rid) = sum(abs(thisReact_A.sel_L)==1);
%         RipplesTable_p.n_mfs_L(rid) = sum(abs(thisReact_A.sel_L)==2);
%         RipplesTable_p.n_hom_L(rid) = sum(abs(thisReact_A.sel_L)==3);
%         RipplesTable_p.n_het_L(rid) = sum(abs(thisReact_A.sel_L)==4);
% 
%         RipplesTable_p.n_sf_R(rid) = sum(abs(thisReact_A.sel_R)==1);
%         RipplesTable_p.n_mfs_R(rid) = sum(abs(thisReact_A.sel_R)==2);
%         RipplesTable_p.n_hom_R(rid) = sum(abs(thisReact_A.sel_R)==3);
%         RipplesTable_p.n_het_R(rid) = sum(abs(thisReact_A.sel_R)==4);
% 
%         RipplesTable_p.n_sf_C(rid) = sum(abs(thisReact_A.sel_C)==1);
%         RipplesTable_p.n_mfs_C(rid) = sum(abs(thisReact_A.sel_C)==2);
%         RipplesTable_p.n_hom_C(rid) = sum(abs(thisReact_A.sel_C)==3);
%         RipplesTable_p.n_het_C(rid) = sum(abs(thisReact_A.sel_C)==4);
% 
%         temp = [sum(thisReact_B.sel_L<=-0.1),sum(thisReact_B.sel_L>=0.1)]; [m,t1] = max(temp); [n,t2] = min(temp);
%         RipplesTable_p.n_fields_L(rid) = (-1)^t1 * m;
%         RipplesTable_p.n_afields_L(rid) = (-1)^t2 * n;
% 
%         
%         temp = [sum(thisReact_B.sel_R<=-0.1),sum(thisReact_B.sel_R>=0.1)]; [m,t1] = max(temp); [n,t2] = min(temp);
%         RipplesTable_p.n_fields_R(rid)= (-1)^t1 * m;
%         RipplesTable_p.n_afields_R(rid) = (-1)^t2 * n;
% 
%         temp = [sum(thisReact_B.sel_C<=-0.1),sum(thisReact_B.sel_C>=0.1)]; [m,t1] = max(temp); [n,t2] = min(temp);
%         RipplesTable_p.n_fields_C(rid)= (-1)^t1 * m;
%         RipplesTable_p.n_afields_C(rid) = (-1)^t2 * n;
% 
% 
% 
%     ensi = [RDI_L.act_mean,RDI_R.act_mean,RDI_C.act_mean];
%     [m,t1] = max(abs(ensi));
% 
%     RipplesTable_p.mRDI_M_UV(rid) = ensi(t1);

catch

    thisUnits=table;
    for u=1:size(thisReact_B,1)
        thisUnits =[thisUnits;UnitsTable(find(cellfun(Params.cellfind(thisReact_B.UnitID(u)),UnitsTable.ID)),:)];
    end


    if ~isempty(thisUnits)
            RipplesTable_p.mRDI_L_UV(rid) = mean(thisUnits.RDI_LScene(abs(thisUnits.RDI_LScene)>=0.1));
        RipplesTable_p.mRDI_R_UV(rid) = mean(thisUnits.RDI_RScene(abs(thisUnits.RDI_RScene)>=0.1));
        RipplesTable_p.mRDI_C_UV(rid) = mean(thisUnits.RDI_LR(abs(thisUnits.RDI_LR)>=0.1));

        if isnan(RipplesTable_p.mRDI_L_UV(rid)), RipplesTable_p.mRDI_L_UV(rid) = nanmean(thisUnits.RDI_LScene); end
        if isnan(RipplesTable_p.mRDI_R_UV(rid)), RipplesTable_p.mRDI_R_UV(rid) = nanmean(thisUnits.RDI_RScene); end
        if isnan(RipplesTable_p.mRDI_C_UV(rid)), RipplesTable_p.mRDI_C_UV(rid) = nanmean(thisUnits.RDI_LR); end

               RipplesTable_p.snr_RDI_L_UV(rid) = cal_snr(thisUnits.RDI_LScene);
       RipplesTable_p.snr_RDI_R_UV(rid) = cal_snr(thisUnits.RDI_RScene); 
       RipplesTable_p.snr_RDI_C_UV(rid) = cal_snr(thisUnits.RDI_LR);
    else
               RipplesTable_p.snr_RDI_L_UV(rid) = nan;
       RipplesTable_p.snr_RDI_R_UV(rid) = nan; 
       RipplesTable_p.snr_RDI_C_UV(rid) = nan;

           RipplesTable_p.WilRDI_L_UV(rid) = nan;
    RipplesTable_p.MRDI_L_UV(rid) = nan;
     RipplesTable_p.pBinom_L_UV(rid) = nan;

           RipplesTable_p.WilRDI_R_UV(rid) = nan;
    RipplesTable_p.MRDI_R_UV(rid) = nan;
     RipplesTable_p.pBinom_R_UV(rid) = nan;

         RipplesTable_p.WilRDI_C_UV(rid) = nan;
    RipplesTable_p.MRDI_C_UV(rid) = nan;
     RipplesTable_p.pBinom_C_UV(rid) = nan;

    end

    RipplesTable_p.MRDI_L_UV(rid) = nan;
RipplesTable_p.MRDI_R_UV(rid) = nan;
RipplesTable_p.MRDI_C_UV(rid) = nan;
    end

end
% RipplesTable_p = RipplesTable_p(RipplesTable_p.nFields>=5,:);
suff = '_forAnalysis';
writetable(RipplesTable_p,[ROOT.Save '\RipplesTable_' thisRegion0 suff '.xlsx'],'writemode','replacefile')
% delete(gcp('nocreate'));



%%
% RipplesTable_p.snr_RDI_L_UV(RipplesTable_p.snr_RDI_L_UV==20) =25;
% RipplesTable_p.snr_RDI_L_UV(RipplesTable_p.snr_RDI_L_UV==-20) =-25;
%%



function [x1,x2] = rand_samp(thisPool_A,thisPool_B,thisUnits,rdi)
samp = datasample(thisPool_A,size(thisUnits,1),1);
if ~isempty(samp)
    samp2=table;
    for j=1:size(thisPool_B,1)
        if max(strcmp(thisPool_B.ID{j}(1:12),samp.ID)) & abs(thisPool_B.(rdi)(j))>0.1
            samp2=[samp2;thisPool_B(j,:)];
        end
    end
    x1 = nanmean(samp2.(rdi)); x2 = nanmedian(samp2.(rdi));
else
    x1 = nan; x2=nan;
end
end
