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
%
SessionTable  = readtable([ROOT.Save '\SessionCountTable_' thisRegion '_forAnalysis_RDI.xlsx']);

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
    Trials.choice_correct = mod(Trials.context,2);
    Trials.choice_correct(Trials.choice_correct==0)=2;
    thisSripples = RipplesTable(RipplesTable.rat==SessionTable.rat(sid) & RipplesTable.session==SessionTable.session(sid),:);

    cid = Trials.correctness==1;

    SessionTable.Wrong_ITI(sid) = sum(Trials.duration_ITI(find(~cid)));
    SessionTable.Correct_ITI(sid) = sum(Trials.duration_ITI(find(cid)));

    cid(1)=1;
    SessionTable.Bef_Wrong_ITI(sid) = sum(Trials.duration_ITI(find(~cid)-1));

    cid(1)=0;
    SessionTable.Bef_Correct_ITI(sid) = sum(Trials.duration_ITI(find(cid)-1));

    cid = Trials.correctness==1;
SessionTable.All_Rips(sid) = size(thisSripples.tr_temp,1);

SessionTable.Correct_Rips(sid) = sum(ismember(thisSripples.tr_temp,find(cid)));
SessionTable.Wrong_Rips(sid) = sum(ismember(thisSripples.tr_temp,find(~cid)));

% SessionTable.Bef_Correct_Rips(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)));
% SessionTable.Bef_Wrong_Rips(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)));
% 

SessionTable.Correct_Replays(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) & thisSripples.Decoding_All);
SessionTable.Wrong_Replays(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))& thisSripples.Decoding_All);

% SessionTable.Bef_Correct_Replays(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) & thisSripples.Decoding_All);
% SessionTable.Bef_Wrong_Replays(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) & thisSripples.Decoding_All);


SessionTable.Correct_SceneB(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) & thisSripples.SceneBias);
SessionTable.Wrong_SceneB(sid) = sum(ismember(thisSripples.tr_temp,find(~cid)) & thisSripples.SceneBias);

% SessionTable.Bef_Correct_SceneB(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) & thisSripples.SceneBias);
% SessionTable.Bef_Wrong_SceneB(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) & thisSripples.SceneBias);
% 

SessionTable.Correct_ChoiceB(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) & thisSripples.ChoiceBias);
SessionTable.Wrong_ChoiceB(sid) = sum(ismember(thisSripples.tr_temp,find(~cid)) & thisSripples.ChoiceBias);

% SessionTable.Bef_Correct_ChoiceB(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) & thisSripples.ChoiceBias);
% SessionTable.Bef_Wrong_ChoiceB(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) & thisSripples.ChoiceBias);
% 

SessionTable.Correct_SReplays(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) & thisSripples.Decoding_All & thisSripples.SceneBias);
SessionTable.Wrong_SReplays(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))& thisSripples.Decoding_All & thisSripples.SceneBias);

% SessionTable.Bef_Correct_SReplays(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) & thisSripples.Decoding_All & thisSripples.SceneBias);
% SessionTable.Bef_Wrong_SReplays(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) & thisSripples.Decoding_All & thisSripples.SceneBias);
% 

SessionTable.Correct_CReplays(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) & thisSripples.Decoding_All & thisSripples.ChoiceBias);
SessionTable.Wrong_CReplays(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))& thisSripples.Decoding_All & thisSripples.ChoiceBias);

% SessionTable.Bef_Correct_CReplays(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) & thisSripples.Decoding_All & thisSripples.ChoiceBias);
% SessionTable.Bef_Wrong_CReplays(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) & thisSripples.Decoding_All & thisSripples.ChoiceBias);

Trials.context(length(cid)+1)=nan;
SessionTable.Correct_SMatch(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) &...
    (thisSripples.SceneBias==Trials.context(thisSripples.tr_temp)));
SessionTable.Wrong_SMatch(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))&...
    (thisSripples.SceneBias==Trials.context(thisSripples.tr_temp)));
% SessionTable.Bef_Correct_SMatch(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) &...
%     (thisSripples.SceneBias==Trials.context(thisSripples.tr_temp+1)));
% SessionTable.Bef_Wrong_SMatch(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) &...
%     (thisSripples.SceneBias==Trials.context(thisSripples.tr_temp+1)));

SessionTable.Correct_CMatch(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) &...
    (thisSripples.ChoiceBias==Trials.choice_side(thisSripples.tr_temp)));
SessionTable.Wrong_CMatch(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))&...
    (thisSripples.ChoiceBias==Trials.choice_side(thisSripples.tr_temp)));
% SessionTable.Bef_Correct_CMatch(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) &...
%     (thisSripples.ChoiceBias==Trials.choice_side(thisSripples.tr_temp+1)));
% SessionTable.Bef_Wrong_CMatch(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) &...
%     (thisSripples.ChoiceBias==Trials.choice_side(thisSripples.tr_temp+1)));

SessionTable.Correct_CMatch_correct(sid) = sum(ismember(thisSripples.tr_temp,find(cid)) &...
   (thisSripples.ChoiceBias==Trials.choice_correct(thisSripples.tr_temp)));
SessionTable.Wrong_CMatch_correct(sid) = sum(ismember(thisSripples.tr_temp,find(~cid))&...
    (thisSripples.ChoiceBias==Trials.choice_correct(thisSripples.tr_temp)));
% SessionTable.Bef_Correct_CMatch_correct(sid) = sum(ismember(thisSripples.tr_temp+1,find(cid)) &...
%     (thisSripples.ChoiceBias==Trials.choice_correct(thisSripples.tr_temp+1)));
% SessionTable.Bef_Wrong_CMatch_correct(sid) = sum(ismember(thisSripples.tr_temp+1,find(~cid)) &...
%     (thisSripples.ChoiceBias==Trials.choice_correct(thisSripples.tr_temp+1)));
end

% writetable(SessionTable,[ROOT.Save '\SessionCountTable_' thisRegion '_forAnalysis_RDI.xlsx']);


%%
S = SessionTable;

%%
figure;

subplot(2,3,1)
bar4(S.Correct_Rips./S.Bef_Correct_ITI,...
S.Bef_Wrong_Rips./S.Bef_Wrong_ITI,...
S.Correct_Rips./S.Correct_ITI,...
S.Wrong_Rips./S.Wrong_ITI,'Ripple rate (Hz)')

subplot(2,3,2)
bar4(S.Bef_Correct_Replays./S.Bef_Correct_Rips,...
S.Bef_Wrong_Replays./S.Bef_Wrong_Rips,...
S.Correct_Replays./S.Correct_Rips,...
S.Wrong_Replays./S.Wrong_Rips, 'Replay rate')

subplot(2,3,3)
bar4(S.Bef_Correct_SceneB./S.Bef_Correct_Rips,...
S.Bef_Wrong_SceneB./S.Bef_Wrong_Rips,...
S.Correct_SceneB./S.Correct_Rips,...
S.Wrong_SceneB./S.Wrong_Rips, 'Scene-selectivite ripple rate')

subplot(2,3,4)
bar4(S.Bef_Correct_ChoiceB./S.Bef_Correct_Rips,...
S.Bef_Wrong_ChoiceB./S.Bef_Wrong_Rips,...
S.Correct_ChoiceB./S.Correct_Rips,...
S.Wrong_ChoiceB./S.Wrong_Rips, 'Choice-selectivite ripple rate')

%%
figure;
inters = sum(S.Wrong_SReplays+S.Wrong_CReplays);
xw = [sum([S.Wrong_SceneB;S.Wrong_ChoiceB])-inters,inters,sum(S.Wrong_Replays)-inters,...
    sum(S.Wrong_Rips)-sum([S.Wrong_ChoiceB;S.Wrong_SceneB;S.Wrong_Replays])+inters];
inters = sum(S.Correct_SReplays+S.Correct_CReplays);
xc = [sum([S.Correct_ChoiceB;S.Correct_SceneB])-inters, inters,sum(S.Correct_Replays)-inters,...
    sum(S.Correct_Rips)-sum([S.Correct_ChoiceB;S.Correct_SceneB;S.Correct_Replays])+inters];

 t1 = (S.Correct_Rips)-([S.Correct_ChoiceB+S.Correct_SceneB+S.Correct_Replays])+(S.Correct_SReplays+S.Correct_CReplays);

t2 = (S.Wrong_Rips)-([S.Wrong_ChoiceB+S.Wrong_SceneB+S.Wrong_Replays])+(S.Wrong_SReplays+S.Wrong_CReplays);

bar(flip([xc/sum(xc);xw/sum(xw)],2),'stacked')
[p,h]=ttest2(t1./S.Correct_Rips, t2./S.Wrong_Rips);



figure;
xw = [sum([S.Wrong_SMatch]),sum(S.Wrong_SceneB)-sum([S.Wrong_SMatch])];
xc = [sum([S.Correct_SMatch]),sum(S.Correct_SceneB)-sum([S.Correct_SMatch])];

bar(flip([xc/sum(xc);xw/sum(xw)],2),'stacked')

[p,h]=ttest(S.Wrong_SMatch./S.Wrong_SceneB,S.Correct_SMatch./S.Correct_SceneB);


figure;
xw = [sum([S.Wrong_CMatch]),sum([S.Wrong_CMatch_correct])];
xc = [sum([S.Correct_CMatch]),sum(S.Correct_ChoiceB)-sum([S.Correct_CMatch])];

bar(flip([xc/sum(xc);xw/sum(xw)],2),'stacked')

[p,h]=ranksum(S.Wrong_CMatch./S.Wrong_ChoiceB,S.Correct_CMatch./S.Correct_ChoiceB);

%% compare correct vs. wrong ripples
RC = RipplesTable(RipplesTable.correctness==1,:);
RW = RipplesTable(RipplesTable.correctness==2,:);
%%
figure('position',[680,462,459,516]);
bar2(RC.RippleDuration,RW.RippleDuration,'Ripple duration (s)')

figure('position',[680,462,459,516]);
bar2(RC.spike./(RC.RippleDuration),RW.spike./(RW.RippleDuration),'spike rate (spike / sec)')

figure('position',[680,462,459,516]);
bar2(RC.ensemble,RW.ensemble,'ensemble')

figure('position',[680,462,459,516]);
bar2(RC.MeanFreq,RW.MeanFreq,'Frequency (Hz)')

%%
figure('position',[680,462,459,516]);
bar2((RC.abs_mRDI_C),abs(RW.abs_mRDI_C),'Choice selectivity (abs mean)')

figure('position',[680,462,459,516]);
bar2((RC.mRDI_C),(RW.mRDI_C),'Choice selectivity (mean)')
%%

figure('position',[680,462,459,516]);
RC = RC(RC.nnRDI_C~=0,:);
bar2(abs(((RC.npRDI_C+RC.nnRDI_C)./ RC.nRDIsMax)),abs(((RW.npRDI_C+RW.nnRDI_C)./ RW.nRDIsMax)),'Choice selective cell rate')
%%
figure;

subplot(1,2,1)
bar2(S.Correct_Rips./S.Correct_ITI,...
S.Wrong_Rips./S.Wrong_ITI,'Ripple rate (Hz)')

subplot(1,2,2)
bar2((S.Correct_Replays)./S.Correct_Rips,...
(S.Wrong_Replays)./S.Wrong_Rips, 'Replay rate')


figure;
subplot(2,3,1)
bar2((S.Correct_SceneB+S.Correct_ChoiceB)./S.Correct_Rips,...
(S.Wrong_SceneB+S.Wrong_ChoiceB)./S.Wrong_Rips, 'Scene/Choice selectivite ripple rate')
ylim([0 .3])

subplot(2,3,2)
bar2((S.Correct_SceneB)./S.Correct_Rips,...
(S.Wrong_SceneB)./S.Wrong_Rips, 'Scene-selectivite ripple rate')
ylim([0 .3])

subplot(2,3,3)
bar2((S.Correct_ChoiceB)./S.Correct_Rips,...
(S.Wrong_ChoiceB)./S.Wrong_Rips, 'Choice-selectivite ripple rate')
ylim([0 .3])


subplot(2,3,4)
bar2((S.Correct_SceneB+S.Correct_ChoiceB-S.Correct_SReplays-S.Correct_CReplays)./S.Correct_Rips,...
(S.Wrong_SceneB+S.Wrong_ChoiceB-S.Wrong_SReplays-S.Wrong_CReplays)./S.Wrong_Rips, 'Scene/Choice selectivite ripple rate')
ylim([0 .3])

subplot(2,3,5)
bar2((S.Correct_SceneB-S.Correct_SReplays)./S.Correct_Rips,...
(S.Wrong_SceneB-S.Wrong_SReplays)./S.Wrong_Rips, 'Scene-selectivite ripple rate')
ylim([0 .3])

subplot(2,3,6)
bar2((S.Correct_ChoiceB-S.Correct_CReplays)./S.Correct_Rips,...
(S.Wrong_ChoiceB-S.Wrong_CReplays)./S.Wrong_Rips, 'Choice-selectivite ripple rate')
ylim([0 .3])


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

%%
function bar2(x1,x2,ylab)
% x1=x1(x1~=0); x2=x2(x2~=0); 
X = [nanmean(x1),nanmean(x2)];
b = bar(X);
 b.FaceColor = 'flat';
 hold on
 err = [nanstd(x1)/sqrt(sum(~isnan(x1))) nanstd(x2)/sqrt(sum(~isnan(x2)))];

 er = errorbar([1 2], X, err,-err);
 er.Color=[0 0 0];
%  er.LineStyle = 'none';
 er.LineWidth = 2;
%  b.CData(3,:)= [1 0 0];
%   b.CData(4,:)= [1 0 0];

% ylim([0 max(X)*1.1])
xticklabels({'Correct','Wrong'})
ylabel(ylab)

[p,~] = ranksum(x1,x2);
text(1.2,max(X)*1.05,['p=' jjnum2str(p,3)],'fontsize',12,'fontweight','b')


set(gca,'fontsize',12,'fontweight','b')
end
