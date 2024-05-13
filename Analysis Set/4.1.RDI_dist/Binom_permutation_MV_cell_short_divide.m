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

RegionList = {'SUB','CA1'};
for reg=1:2
    for crit_si=.05:.05:.4
thisR = RegionList{reg};
% 

thisRegion0 = thisR;
thisRegion = thisR;
thisRegion2 = [thisR '_field'];

% thisRegion0 = 'CA1';
% thisRegion = 'CA1';
% thisRegion2 = 'CA1_field';

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
UnitsTableA = UnitsTable_A;

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

        thisUnitsA=table;
    for u=1:size(thisReact_A,1)
        thisUnitsA =[thisUnitsA;UnitsTableA(find(cellfun(Params.cellfind(thisReact_A.UnitID(u)),UnitsTableA.ID)),:)];
    end

%         RipplesTable_p.mRDI_L_UV(rid) = RDI_L.act_mean;
%         RipplesTable_p.mRDI_R_UV(rid) = RDI_R.act_mean;
%         RipplesTable_p.mRDI_C_UV(rid) = RDI_C.act_mean;



%%
[m,b,b1,b2,thisUnits_p,thisUnits_n] = FiltUnits(thisUnits,thisUnitsA,UnitsTableA,thisReact_A, UnitsTable_B,'RDI_LScene',crit_si);
p1=0;p2=0;
% q1 = sum(thisUnitsA.RDI_LScene>=0.1)/sum(abs(thisUnitsA.RDI_LScene)>=0.1); q2 = sum(thisUnitsA.RDI_LScene<=-0.1)/sum(abs(thisUnitsA.RDI_LScene)>=0.1);
q1 = sum(thisUnits_p.RDI_LScene>=crit_si)/sum(abs(thisUnits_p.RDI_LScene)>=crit_si); if(sum(abs(thisUnits_p.RDI_LScene)>=crit_si))<5, q1=nan;end
q2 = sum(thisUnits_n.RDI_LScene<=-crit_si)/sum(abs(thisUnits_n.RDI_LScene)>=crit_si); if(sum(abs(thisUnits_n.RDI_LScene)>=crit_si))<5, q2=nan;end
RipplesTable_p.pRatio_L_UV(rid) = max([q1-p1,q2-p2]);
RipplesTable_p.pBinomDev_L_UV(rid) = min([b1,b2]);
if RipplesTable_p.pRatio_L_UV(rid)>=0.8, RipplesTable_p.pBinomDev_L_UV(rid) = min([b,b1,b2]); end
RipplesTable_p.nRDI_L_max(rid) = sum(abs(thisUnitsA.RDI_LScene)>=crit_si);

[m,b,b1,b2,thisUnits_p,thisUnits_n] = FiltUnits(thisUnits,thisUnitsA,UnitsTableA,thisReact_A, UnitsTable_B,'RDI_RScene',crit_si);
p1=0;p2=0;
% q1 = sum(thisUnitsA.RDI_RScene>=0.1)/sum(abs(thisUnitsA.RDI_RScene)>=0.1); q2 = sum(thisUnitsA.RDI_RScene<=-0.1)/sum(abs(thisUnitsA.RDI_RScene)>=0.1);
q1 = sum(thisUnits_p.RDI_RScene>=crit_si)/sum(abs(thisUnits_p.RDI_RScene)>=crit_si); if(sum(abs(thisUnits_p.RDI_RScene)>=crit_si))<5, q1=nan;end
q2 = sum(thisUnits_n.RDI_RScene<=-crit_si)/sum(abs(thisUnits_n.RDI_RScene)>=crit_si); if(sum(abs(thisUnits_n.RDI_RScene)>=crit_si))<5, q2=nan;end
RipplesTable_p.pRatio_R_UV(rid) = max([q1-p1,q2-p2]);
RipplesTable_p.pBinomDev_R_UV(rid) = min([b1,b2]);
if RipplesTable_p.pRatio_R_UV(rid)>=0.8, RipplesTable_p.pBinomDev_R_UV(rid) = min([b,b1,b2]); end
RipplesTable_p.nRDI_R_max(rid) = sum(abs(thisUnitsA.RDI_RScene)>=crit_si);

[m,b,b1,b2,thisUnits_p,thisUnits_n] = FiltUnits(thisUnits,thisUnitsA,UnitsTableA,thisReact_A, UnitsTable_B,'RDI_LR',crit_si);
p1=0;p2=0;
% q1 = sum(thisUnitsA.RDI_LR>=0.1)/sum(abs(thisUnitsA.RDI_LR)>=0.1); q2 = sum(thisUnitsA.RDI_LR<=-0.1)/sum(abs(thisUnitsA.RDI_LR)>=0.1);
q1 = sum(thisUnits_p.RDI_LR>=crit_si)/sum(abs(thisUnits_p.RDI_LR)>=crit_si); if(sum(abs(thisUnits_p.RDI_LR)>=crit_si))<5, q1=nan;end
q2 = sum(thisUnits_n.RDI_LR<=-crit_si)/sum(abs(thisUnits_n.RDI_LR)>=crit_si); if(sum(abs(thisUnits_n.RDI_LR)>=crit_si))<5, q2=nan;end
RipplesTable_p.pRatio_C_UV(rid) = max([q1-p1,q2-p2]);
RipplesTable_p.pBinomDev_C_UV(rid) = min([b1,b2]);
if RipplesTable_p.pRatio_C_UV(rid)>=0.8, RipplesTable_p.pBinomDev_C_UV(rid) = min([b,b1,b2]); end
RipplesTable_p.nRDI_C_max(rid) = sum(abs(thisUnitsA.RDI_LR)>=crit_si);
%%
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
        RipplesTable_p.mRDI_L_UV(rid) = mean(thisUnits.RDI_LScene(abs(thisUnits.RDI_LScene)>=crit_si));
        RipplesTable_p.mRDI_R_UV(rid) = mean(thisUnits.RDI_RScene(abs(thisUnits.RDI_RScene)>=crit_si));
        RipplesTable_p.mRDI_C_UV(rid) = mean(thisUnits.RDI_LR(abs(thisUnits.RDI_LR)>=crit_si));

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
% suff = '_forAnalysis_240310';
% writetable(RipplesTable_p,[ROOT.Save '\RipplesTable_' thisRegion0 suff '.xlsx'],'writemode','replacefile')
% delete(gcp('nocreate'));
RipplesTable_SI.(thisRegion0).(['SI_' jmnum2str(crit_si*100,3)]) = RipplesTable_p;
    end
end

save([ROOT.Processed '\Manuscript\ripple_SI_4.mat'],'RipplesTable_SI');
%%
% RipplesTable_p.snr_RDI_L_UV(RipplesTable_p.snr_RDI_L_UV==20) =25;
% RipplesTable_p.snr_RDI_L_UV(RipplesTable_p.snr_RDI_L_UV==-20) =-25;
%%


%%
function [m,b,b1,b2,thisUnits_p,thisUnits_n] = FiltUnits(thisUnits,thisUnitsA,UnitsTableA,thisReact_A, UnitsTable_B,var,crit_si)

    thisUnits_p=table; thisUnits_n=table;
    for u=1:size(thisReact_A,1)
        temp =UnitsTable_B(find(strncmp(thisReact_A.UnitID(u),UnitsTable_B.ID,12)),:);
        if ~isempty(temp)
            if ~isnan(max(temp.(var)))
        [~,m] = max((temp.(var)));
        thisUnits_p =[thisUnits_p;temp(m,:)];

         [~,m] = min((temp.(var)));
        thisUnits_n =[thisUnits_n;temp(m,:)];
            else
                thisUnits_p =[thisUnits_p;temp(1,:)];
                thisUnits_n =[thisUnits_n;temp(1,:)];
            end
        end

    end

    m = median(UnitsTable_B.(var)(abs(UnitsTable_B.(var))>=crit_si));
p1 = sum(UnitsTable_B.(var)>=crit_si)/sum(abs(UnitsTable_B.(var))>=crit_si); p2 = sum(UnitsTable_B.(var)<=-crit_si)/sum(abs(UnitsTable_B.(var))>=crit_si);
p00 = sum(UnitsTableA.(var)>=crit_si)/sum(abs(UnitsTableA.(var))>=crit_si); p01 = sum(UnitsTableA.(var)<=-crit_si)/sum(abs(UnitsTableA.(var))>=crit_si);
 

if sum(abs(thisUnits.(var))>=crit_si)<1
b=nan; b1 = nan; b2=nan;
else
    b00=myBinomTest(sum(thisUnitsA.(var)>=crit_si),sum(abs(thisUnitsA.(var))>=crit_si),p2,'one');
     b01=myBinomTest(sum(thisUnitsA.(var)<=-crit_si),sum(abs(thisUnitsA.(var))>=crit_si),p1,'one');

     b = min([b00,b01]);
%      b=1;
b1=myBinomTest(sum(thisUnits_p.(var)>=crit_si),sum(abs(thisUnits_p.(var))>=crit_si),p1,'one');
b2= myBinomTest(sum(thisUnits_n.(var)<=-crit_si),sum(abs(thisUnits_n.(var))>=crit_si),p2,'one');
% b3= myBinomTest(sum(thisUnits_n.(var)<=-crit_si),sum(abs(thisUnits_n.(var))>=crit_si),p1,'one');
    end
end