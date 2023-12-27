Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip1 = [ROOT.Save '\ripples_mat\R1'];
ROOT.Rip2 = [ROOT.Save '\ripples_mat\R2'];
ROOT.Rip3 = [ROOT.Save '\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Save '\ripples_mat\R4'];
ROOT.Rip5 = [ROOT.Save '\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Save '\Manuscript figures\R2\R0_fig'];
ROOT.Units = [ROOT.Save '\units_mat\U2'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];

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

SessionRips=table;
sid=1;
%%
thisRegion0 = 'SUB';

RA = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis.xlsx']);
R0 = readtable([ROOT.Rip0 '\RipplesList_' thisRegion0 '.xlsx']);
R2 = readtable([ROOT.Rip2 '\RipplesTable_Behav_' thisRegion0 '_speed_filtered.xlsx']);
R3 = readtable([ROOT.Rip3 '\RipplesTable_Behav_Overlap' thisRegion0 '.xlsx']);
R4 = readtable([ROOT.Rip3 '\RipplesTable_' thisRegion0 '_forAnalysis.xlsx']);
% RipplesTable = readtable([ROOT.Rip '\RipplesTable_' thisRegion0  '_forAnalysis.xlsx']);

TT_table = readtable([ROOT.Info '\TT_table.xlsx']);


for s=1:size(SessionList,1)
    % for sid=32:32
    if ~SessionList.include(s) | ~strcmp(SessionList.experimenter(s),'LSM'), continue; end
    thisRIDn = SessionList.rat(s); thisRID=jmnum2str(thisRIDn,3);
    thisSIDn = SessionList.session(s); thisSID=jmnum2str(thisSIDn,2);
    thisRSID = [thisRID '-' thisSID];

    BehavTable = readtable([ROOT.Behav '\' thisRSID '.xlsx']);


    theseR0 = R0(R0.rat==thisRIDn & R0.session==thisSIDn,:);
    theseR2 = R2(R2.rat==thisRIDn & R2.session==thisSIDn,:);
    theseRA = RA(RA.rat==thisRIDn & RA.session==thisSIDn,:);
    
    SessionRips.session{sid} = thisRSID;
    SessionRips.region{sid} = thisRegion0;
    SessionRips.time_sum(sid) = sum(BehavTable.duration_ITI(BehavTable.correctness==1));

    SessionRips.Rips_LFP(sid) = size(theseR0,1);
    SessionRips.Rips_Overlap(sid) = sum(theseR0.Overlap~=0);
    SessionRips.Rips_Ensemble(sid) = sum(theseR0.ensemble>2);
    SessionRips.Rips_O_E(sid) = sum(theseR0.Overlap~=0 & theseR0.ensemble>2);
    SessionRips.Rips_O_E_Behav(sid) = sum(theseR2.Overlap~=0 & theseR2.ensemble>2);
    SessionRips.Rips_O_E_B_Pc(sid) = sum(theseRA.Overlap~=0);

    sid=sid+1;
end
%%
writetable(SessionRips,[ROOT.Save '\SessionRips_pre.xlsx'],'WriteMode','replacefile')

%  SessionRips.RippleRate = SessionRips.num_rips./SessionRips.time_sum;
%%
SessionRips_sub = SessionRips;
SessionRips_sub(isnan(SessionRips_sub.RippleRate),:)=[];
SessionRips_ca1(isnan(SessionRips_ca1.RippleRate),:)=[];


 s0 = SessionRips_sub.num_rips;
  s1 = SessionRips_ca1.num_rips;

   nanmean(s0)
  nanstd(s0)

 nanmean(s1)
  nanstd(s1)

 [h,p]= ttest2(s1,s0);