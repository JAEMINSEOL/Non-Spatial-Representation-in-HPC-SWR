Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Save '\ripples_mat\R4_10ms'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R11_AI_hot'];
ROOT.Units = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];

dir = '';
ROOT.Fig = [ROOT.Fig3];
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 1;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.05*Params.Fs;
dur = 0.4*Params.Fs;
thisFRMapSCALE=2;
Params.tbinDuration = 0.005;

temp = RipplesTable.rat+RipplesTable.session*10^(-2);
[n,t] = unique(temp);
SessionTable = table;
SessionTable.rat = RipplesTable.rat(t);
SessionTable.session = RipplesTable.session(t);

for rid = 1:size(RipplesTable,1)
    RipplesTable.tr_temp(rid) = str2double(RipplesTable.trial{rid}(end-2:end));
end

%%
for sid = 1:size(SessionTable,1)
    Trials = readtable([ROOT.Behav '\' jmnum2str(SessionTable.rat(sid),3) '-' jmnum2str(SessionTable.session(sid),2) '.xlsx']);
    thisSripples = RipplesTable(RipplesTable.rat==SessionTable.rat(sid) & RipplesTable.session==SessionTable.session(sid),:);

    cid = Trials.correctness==1;

    SessionTable.Aft_Wrong_ITI(sid) = sum(Trials.duration_ITI(find(~cid)));
    SessionTable.Aft_Correct_ITI(sid) = sum(Trials.duration_ITI(find(cid)));

    cid(1)=1;
    SessionTable.Bef_Wrong_ITI(sid) = sum(Trials.duration_ITI(find(~cid)-1));

    cid(1)=0;
    SessionTable.Bef_Correct_ITI(sid) = sum(Trials.duration_ITI(find(cid)-1));

    cid = Trials.correctness==1;
SessionTable.All_Rips(sid) = size(thisSripples.tr_temp,1);

SessionTable.Aft_Correct_Rips(sid) = sum(ismember(thisSripples.tr_temp,find(cid)));
SessionTable.Aft_Wrong_Rips(sid) = sum(ismember(thisSripples.tr_temp,find(~cid)));

SessionTable.Bef_Correct_Rips(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)));
SessionTable.Bef_Wrong_Rips(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)));


SessionTable.Aft_Correct_Replays(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) & thisSripples.Decoding_All);
SessionTable.Aft_Wrong_Replays(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))& thisSripples.Decoding_All);

SessionTable.Bef_Correct_Replays(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) & thisSripples.Decoding_All);
SessionTable.Bef_Wrong_Replays(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) & thisSripples.Decoding_All);


SessionTable.Aft_Correct_SceneB(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) & thisSripples.SceneBias);
SessionTable.Aft_Wrong_SceneB(sid) = sum(ismember(thisSripples.tr_temp,find(~cid)) & thisSripples.SceneBias);

SessionTable.Bef_Correct_SceneB(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) & thisSripples.SceneBias);
SessionTable.Bef_Wrong_SceneB(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) & thisSripples.SceneBias);


SessionTable.Aft_Correct_ChoiceB(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) & thisSripples.ChoiceBias);
SessionTable.Aft_Wrong_ChoiceB(sid) = sum(ismember(thisSripples.tr_temp,find(~cid)) & thisSripples.ChoiceBias);

SessionTable.Bef_Correct_ChoiceB(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) & thisSripples.ChoiceBias);
SessionTable.Bef_Wrong_ChoiceB(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) & thisSripples.ChoiceBias);


SessionTable.Aft_Correct_SReplays(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) & thisSripples.Decoding_All & thisSripples.SceneBias);
SessionTable.Aft_Wrong_SReplays(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))& thisSripples.Decoding_All & thisSripples.SceneBias);

SessionTable.Bef_Correct_SReplays(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) & thisSripples.Decoding_All & thisSripples.SceneBias);
SessionTable.Bef_Wrong_SReplays(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) & thisSripples.Decoding_All & thisSripples.SceneBias);


SessionTable.Aft_Correct_CReplays(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) & thisSripples.Decoding_All & thisSripples.ChoiceBias);
SessionTable.Aft_Wrong_CReplays(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))& thisSripples.Decoding_All & thisSripples.ChoiceBias);

SessionTable.Bef_Correct_CReplays(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) & thisSripples.Decoding_All & thisSripples.ChoiceBias);
SessionTable.Bef_Wrong_CReplays(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) & thisSripples.Decoding_All & thisSripples.ChoiceBias);

Trials.context(length(cid)+1)=nan;
SessionTable.Aft_Correct_SMatch(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) &...
    thisSripples.Decoding_All & (thisSripples.SceneBias==Trials.context(thisSripples.tr_temp)));
SessionTable.Aft_Wrong_SMatch(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))&...
    thisSripples.Decoding_All & (thisSripples.SceneBias==Trials.context(thisSripples.tr_temp)));
SessionTable.Bef_Correct_SMatch(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) &...
    thisSripples.Decoding_All & (thisSripples.SceneBias==Trials.context(thisSripples.tr_temp+1)));
SessionTable.Bef_Wrong_SMatch(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) &...
    thisSripples.Decoding_All & (thisSripples.SceneBias==Trials.context(thisSripples.tr_temp+1)));

SessionTable.Aft_Correct_CMatch(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) &...
    thisSripples.Decoding_All & (thisSripples.ChoiceBias==Trials.context(thisSripples.tr_temp)));
SessionTable.Aft_Wrong_CMatch(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))&...
    thisSripples.Decoding_All & (thisSripples.ChoiceBias==Trials.context(thisSripples.tr_temp)));
SessionTable.Bef_Correct_CMatch(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) &...
    thisSripples.Decoding_All & (thisSripples.ChoiceBias==Trials.context(thisSripples.tr_temp+1)));
SessionTable.Bef_Wrong_CMatch(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) &...
    thisSripples.Decoding_All & (thisSripples.ChoiceBias==Trials.context(thisSripples.tr_temp+1)));
end

writetable(SessionTable,[ROOT.Save '\SessionCountTable_' thisRegion '_forAnalysis_RDI.xlsx']);



%%
S = SessionTable;
figure;

subplot(2,3,1)
bar4(S.Bef_Correct_Rips./S.Bef_Correct_ITI,...
S.Bef_Wrong_Rips./S.Bef_Wrong_ITI,...
S.Aft_Correct_Rips./S.Aft_Correct_ITI,...
S.Aft_Wrong_Rips./S.Aft_Wrong_ITI,'Ripple rate (Hz)')

subplot(2,3,2)
bar4(S.Bef_Correct_Replays./S.Bef_Correct_Rips,...
S.Bef_Wrong_Replays./S.Bef_Wrong_Rips,...
S.Aft_Correct_Replays./S.Aft_Correct_Rips,...
S.Aft_Wrong_Replays./S.Aft_Wrong_Rips, 'Replay rate')

subplot(2,3,3)
bar4(S.Bef_Correct_SceneB./S.Bef_Correct_Rips,...
S.Bef_Wrong_SceneB./S.Bef_Wrong_Rips,...
S.Aft_Correct_SceneB./S.Aft_Correct_Rips,...
S.Aft_Wrong_SceneB./S.Aft_Wrong_Rips, 'Scene-selectivite ripple rate')

subplot(2,3,4)
bar4(S.Bef_Correct_ChoiceB./S.Bef_Correct_Rips,...
S.Bef_Wrong_ChoiceB./S.Bef_Wrong_Rips,...
S.Aft_Correct_ChoiceB./S.Aft_Correct_Rips,...
S.Aft_Wrong_ChoiceB./S.Aft_Wrong_Rips, 'Choice-selectivite ripple rate')
%%
function bar4(x1,x2,x3,x4,ylab)
x1=x1(x1~=0); x2=x2(x2~=0); x3=x3(x3~=0); x4=x4(x4~=0);
X = [nanmean(x1),nanmean(x2),nanmean(x3),nanmean(x4)];
b = bar(X);
 b.FaceColor = 'flat';
 b.CData(3,:)= [1 0 0];
  b.CData(4,:)= [1 0 0];

ylim([0 max(X)*1.1])
xticklabels({'Before-Correct','Before-Wrong','After-Correct','After-Wrong'})
ylabel(ylab)

[~,p] = ttest2(x1,x2);
text(1.2,max(X)*1.05,['p=' jjnum2str(p,3)],'fontsize',12,'fontweight','b')

[~,p] = ttest2(x3,x4);
text(3.2,max(X)*1.05,['p=' jjnum2str(p,3)],'fontsize',12,'fontweight','b')

set(gca,'fontsize',12,'fontweight','b')
end