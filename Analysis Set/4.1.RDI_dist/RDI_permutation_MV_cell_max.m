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
thisR = RegionList{reg};
% 

thisRegion0 = thisR;
thisRegion = thisR;
thisRegion2 = [thisR '_field'];

% RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis' '.xlsx']);
ReactTable_A = readtable([ROOT.Save '\ReactTable_' thisRegion0 '_' thisRegion '.xlsx']);
ReactTable_B = readtable([ROOT.Save '\ReactTable_' thisRegion0 '_' thisRegion2 '.xlsx']);


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
%%
% UnitsTable_A =UnitsTable;
% sz = size(UnitsTable_A,2);
% sa=29; sb=31;
% col_name_A = string(UnitsTable_A.Properties.VariableNames(sa:sz));
% UnitsTable_B(1,sb:sz+2) = array2table(nan(1,sz-sa+1));
% UnitsTable_B = renamevars(UnitsTable_B, UnitsTable_B.Properties.VariableNames(sb:sz+2),col_name_A);
% for uid=1:size(UnitsTable_B,1)
%     tar = find(strncmp(UnitsTable_A.ID,UnitsTable_B.ID(uid),12));
%     if isempty(tar)
%         UnitsTable_B(uid,sb)={nan};
%         continue;
%     end
%     UnitsTable_B(uid,sb:sz+2) = UnitsTable_A(tar,sa:sz);
%     UnitsTable_B.NumField(uid) = str2double(UnitsTable_B.ID{uid}(end));
% end
% %%
% UnitsTable_B(isnan(UnitsTable_B.NumField),:)=[];
%    writetable(UnitsTable_B,[ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx'],'writemode','replacefile');
UnitsTable = UnitsTable_A;

%%
% RipplesTable_p = RipplesTable(RipplesTable.nPCs>2,:);
RipplesTable_p = RipplesTable;
for rid=1:size(RipplesTable_p,1)
    RipID = RipplesTable_p.ID{rid};
    thisReact_A = ReactTable_A(strcmp(ReactTable_A.RippleID,RipID),:);
    thisReact_A = unique(thisReact_A(:,1:2));

    thisReact_B = ReactTable_B(strcmp(ReactTable_B.RippleID,RipID),:);
    thisReact_B = unique(thisReact_B(:,1:2));

    RipplesTable_p.nPCs(rid) = size(thisReact_A,1);
%     if RipplesTable_p.nPCs(rid)<3, continue; end
    RipplesTable_p.nFields(rid) = size(thisReact_B,1);

    thisUnits=table;
    for u=1:size(thisReact_A,1)
        thisUnits =[thisUnits;UnitsTable(find(cellfun(Params.cellfind(thisReact_A.UnitID(u)),UnitsTable.ID)),:)];
    end


%  RipplesTable_p.nFields(rid) = size(thisUnits,1);
 
    if size(thisUnits,1)<3

        RipplesTable_p.pRDI_L(rid)=nan;
        RipplesTable_p.pRDI_R(rid)=nan;
        RipplesTable_p.pRDI_C(rid)=nan;

        RipplesTable_p.pRDI_L_UV(rid)=nan;
        RipplesTable_p.pRDI_R_UV(rid)=nan;
        RipplesTable_p.pRDI_C_UV(rid)=nan;
        continue;
    end


   
    if RipplesTable_p.nFields(rid)<5
               RipplesTable_p.pRDI_L_UV(rid)=nan;
        RipplesTable_p.pRDI_R_UV(rid)=nan;
        RipplesTable_p.pRDI_C_UV(rid)=nan;
        continue; 
    end

%     thisPool_A = UnitsTable_A(UnitsTable_A.rat==RipplesTable_p.rat(rid) & UnitsTable_A.session==RipplesTable_p.session(rid),:);
%     thisPool_B = UnitsTable_B(UnitsTable_B.rat==RipplesTable_p.rat(rid) & UnitsTable_B.session==RipplesTable_p.session(rid),:);

     thisPool_A = UnitsTable_A;
     thisPool_B = UnitsTable_B;
    %




    thisUnits_Sig_L = thisUnits(abs(thisUnits.RDI_LScene)>0.1,:);
    thisUnits_Sig_R = thisUnits(abs(thisUnits.RDI_RScene)>0.1,:);
    thisUnits_Sig_C = thisUnits(abs(thisUnits.RDI_LR)>0.1,:);

    thisPool_A_Sig_L =  thisPool_A(thisPool_A.Selectivity_LScene>=1 ,:);
    thisPool_A_Sig_R = thisPool_A(thisPool_A.Selectivity_RScene>=1 ,:);
    thisPool_A_Sig_C =  thisPool_A(thisPool_A.Selectivity_LR>=1 ,:);

    %     thisUnits_Sig_L = thisUnits;
    % thisUnits_Sig_R = thisUnits;
    % thisUnits_Sig_C = thisUnits;
    % thisPool_Sig_L = thisPool;
    % thisPool_Sig_R = thisPool;
    % thisPool_Sig_C = thisPool;




    RDI_L.dist=table;
    RDI_R.dist=table;
    RDI_C.dist=table;

    temp1=nan; temp2=nan; temp3=nan; temp4=nan; temp5=nan; temp6=nan;
    temp1_scuv=nan; temp2_scuv=nan; temp3_scuv=nan; temp4_scuv=nan; temp5_scuv=nan; temp6_scuv=nan;

         tic
    try
        parfor i=1:randN
            [x1,x2] = rand_samp(thisPool_A_Sig_L,thisPool_A,thisUnits_Sig_L,'RDI_LScene');
            temp1(i) = x1; temp2(i) = x2;
             [x1,x2] = rand_samp(thisPool_A_Sig_R,thisPool_A,thisUnits_Sig_R,'RDI_RScene');
            temp3(i) = x1; temp4(i) = x2;
            [x1,x2] = rand_samp(thisPool_A_Sig_C,thisPool_A,thisUnits_Sig_C,'RDI_LR');
            temp5(i) = x1; temp6(i) = x2;
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

    %     RDI_L.dist.mean_scuv = temp1_scuv;
    %     RDI_L.dist.median_scuv = temp2_scuv;
    %     RDI_R.dist.mean_scuv = temp3_scuv;
    %     RDI_R.dist.median_scuv = temp4_scuv;
    %     RDI_C.dist.mean_scuv = temp5_scuv;
    %     RDI_C.dist.median_scuv = temp6_scuv;


    RDI_L.act_mean = nanmean(thisUnits_Sig_L.RDI_LScene);
    RDI_R.act_mean = nanmean(thisUnits_Sig_R.RDI_RScene);
    RDI_C.act_mean = nanmean(thisUnits_Sig_C.RDI_LR);
    %     RDI_L.act_mean_scuv = nanmean(thisUnits_UV_SC.RDI_LScene);
    %     RDI_R.act_mean_scuv = nanmean(thisUnits_UV_SC.RDI_RScene);
    %     RDI_C.act_mean_scuv = nanmean(thisUnits_UV_SC.RDI_LR);

    RDI_L.act_median = nanmedian(thisUnits_Sig_L.RDI_LScene);
    RDI_R.act_median = nanmedian(thisUnits_Sig_R.RDI_RScene);
    RDI_C.act_median = nanmedian(thisUnits_Sig_C.RDI_LR);
    %     RDI_L.act_median_scuv = nanmedian(thisUnits_UV_SC.RDI_LScene);
    %     RDI_R.act_median_scuv = nanmedian(thisUnits_UV_SC.RDI_RScene);
    %     RDI_C.act_median_scuv = nanmedian(thisUnits_UV_SC.RDI_LR);

    RDI_L.p_mean = min(sum(RDI_L.act_mean>RDI_L.dist.mean),sum(RDI_L.act_mean<RDI_L.dist.mean))/sum(~isnan(RDI_L.dist.mean));
    RDI_R.p_mean = min(sum(RDI_R.act_mean>RDI_R.dist.mean),sum(RDI_R.act_mean<RDI_R.dist.mean))/sum(~isnan(RDI_R.dist.mean));
    RDI_C.p_mean = min(sum(RDI_C.act_mean>RDI_C.dist.mean),sum(RDI_C.act_mean<RDI_C.dist.mean))/sum(~isnan(RDI_C.dist.mean));

    %     RDI_L.p_mean_scuv = min(sum(RDI_L.act_mean_scuv>RDI_L.dist.mean_scuv),sum(RDI_L.act_mean_scuv<RDI_L.dist.mean_scuv))/sum(~isnan(RDI_L.dist.mean_scuv));
    %     RDI_R.p_mean_scuv = min(sum(RDI_R.act_mean_scuv>RDI_R.dist.mean_scuv),sum(RDI_R.act_mean_scuv<RDI_R.dist.mean_scuv))/sum(~isnan(RDI_R.dist.mean_scuv));
    %     RDI_C.p_mean_scuv = min(sum(RDI_C.act_mean_scuv>RDI_C.dist.mean_scuv),sum(RDI_C.act_mean_scuv<RDI_C.dist.mean_scuv))/sum(~isnan(RDI_C.dist.mean_scuv));

    RDI_L.p_median = min(sum(RDI_L.act_median>RDI_L.dist.median),sum(RDI_L.act_median<RDI_L.dist.median))/...
        sum(~isnan(RDI_L.dist.median));
    RDI_R.p_median = min(sum(RDI_R.act_median>RDI_R.dist.median),sum(RDI_R.act_median<RDI_R.dist.median))/...
        sum(~isnan(RDI_R.dist.median));
    RDI_C.p_median = min(sum(RDI_C.act_median>RDI_C.dist.median),sum(RDI_C.act_median<RDI_C.dist.median))/...
        sum(~isnan(RDI_C.dist.median));

    %     RDI_L.p_median_scuv = min(sum(RDI_L.act_median_scuv>RDI_L.dist.median_scuv),sum(RDI_L.act_median_scuv<RDI_L.dist.median_scuv))/...
    %         sum(~isnan(RDI_L.dist.median_scuv));
    %     RDI_R.p_median_scuv = min(sum(RDI_R.act_median_scuv>RDI_R.dist.median_scuv),sum(RDI_R.act_median_scuv<RDI_R.dist.median_scuv))/...
    %         sum(~isnan(RDI_R.dist.median_scuv));
    %     RDI_C.p_median_scuv = min(sum(RDI_C.act_median_scuv>RDI_C.dist.median_scuv),sum(RDI_C.act_median_scuv<RDI_C.dist.median_scuv))/...
    %         sum(~isnan(RDI_C.dist.median_scuv));

    save([ROOT.Rip ['\' thisRegion0 '-' RipID '_max.mat']],'thisUnits','RDI_L','RDI_R','RDI_C')
    disp([RipID ' is finished!'])


    RipplesTable_p.mRDI_L_max(rid) = RDI_L.act_mean;
    RipplesTable_p.mRDI_R_max(rid) = RDI_R.act_mean;
    RipplesTable_p.mRDI_C_max(rid) = RDI_C.act_mean;

        RipplesTable_p.pRDI_L_max(rid) = RDI_L.p_mean;
    RipplesTable_p.pRDI_R_max(rid) = RDI_R.p_mean;
    RipplesTable_p.pRDI_C_max(rid) = RDI_C.p_mean;



    %         RipplesTable_p.pRDI_L_UV_scuv(r) = RDI_L.p_mean_scuv;
    %     RipplesTable_p.pRDI_R_UV_scuv(r) = RDI_R.p_mean_scuv;
    %     RipplesTable_p.pRDI_C_UV_scuv(r) = RDI_C.p_mean_scuv;
    RipplesTable_p.nRDI_L_max(rid) = size(thisUnits_Sig_L,1);
    RipplesTable_p.nRDI_R_max(rid) = size(thisUnits_Sig_R,1);
    RipplesTable_p.nRDI_C_max(rid) = size(thisUnits_Sig_C,1);


    if size(thisUnits_Sig_L,1)<5, RipplesTable_p.pRDI_L_UV(rid)=nan; end
    if size(thisUnits_Sig_R,1)<5, RipplesTable_p.pRDI_R_UV(rid)=nan; end
    if size(thisUnits_Sig_C,1)<5, RipplesTable_p.pRDI_C_UV(rid)=nan; end

    %     if size(thisUnits_UV_SC,1)<5, RipplesTable_p.pRDI_L_UV_scuv(r)=nan; RipplesTable_p.pRDI_R_UV_scuv(r)=nan; RipplesTable_p.pRDI_C_UV_scuv(r)=nan; end

end
% RipplesTable_p = RipplesTable_p(RipplesTable_p.nFields>=5,:);
suff = '_forAnalysis';
writetable(RipplesTable_p,[ROOT.Save '\RipplesTable_' thisRegion0 suff '.xlsx'],'writemode','replacefile')
delete(gcp('nocreate'));

RipplesTable_p.pRDI_M_UV = nanmin([RipplesTable_p.pRDI_L_UV,RipplesTable_p.pRDI_R_UV,RipplesTable_p.pRDI_C_UV],[],2);
sum(RipplesTable_p.pRDI_M_UV<0.05)
sum(RipplesTable_p.DecodingP_all<0.05)
end
%%


%%



function [x1,x2] = rand_samp(thisPool_A,thisPool_B,thisUnits,rdi)
samp = datasample(thisPool_A,size(thisUnits,1),1);
if ~isempty(samp)

    x1 = nanmean(samp.(rdi)); x2 = nanmedian(samp.(rdi));
else
    x1 = nan; x2=nan;
end
end
%%
