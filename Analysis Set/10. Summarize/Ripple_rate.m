Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Mother '\Processed Data\ripples_mat\R4'];
ROOT.Rip5 = [ROOT.Mother '\Processed Data\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Mother '\Manuscript figures\R2\R0_fig'];
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
thisRegion0 = 'SUB_refCA1';
thisRegion = 'SUB';
thisRegion2 = 'SUB_field';
Sel = 'LScene';
% RipplesTable = readtable([ROOT.Rip '\RipplesTable_' thisRegion '_forAnalysis.xlsx']);
% RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_RDIs_UV_cell_HeteroIn_AllPopul.xlsx']);

RipplesTable = readtable([ROOT.Rip '\RipplesTable_' thisRegion0 '_forAnalysis.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);

TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

for sid=1:size(SessionList,1)
% for sid=32:32
    if ~SessionList.include(sid) | ~strcmp(SessionList.experimenter(sid),'LSM'), continue; end
    thisRIDn = SessionList.rat(sid); thisRID=jmnum2str(thisRIDn,3);
    thisSIDn = SessionList.session(sid); thisSID=jmnum2str(thisSIDn,2);
    thisRSID = [thisRID '-' thisSID];
    load([ROOT.Behav '\' thisRSID '.mat']);
    BehavTable = readtable([ROOT.Behav '\' thisRSID '.xlsx']);
    Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);

    theseRip = RipplesTable(RipplesTable.rat==thisRIDn & RipplesTable.session==thisSIDn,:);

    if isempty(theseRip), continue; end


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

time_sum = sum(BehavTable.duration_ITI(BehavTable.correctness==1))
end