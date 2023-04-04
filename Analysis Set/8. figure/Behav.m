Behav = readtable([ROOT.Behav '/BehavTable.xlsx']);
%%
for t = 1:size(Behav)
    if ismember(str2double(Behav.TrialID{t}(1:3)),[232 234 295 415 561])
    Behav.session(t) = str2double(Behav.TrialID{t}(1:3))*1e2+str2double(Behav.TrialID{t}(5:6));
    else
        Behav.session(t) =nan;
    end
end

%%

CxtPerf = table;
CxtPerf = unique(Behav(:,[15,2]),'rows');

CxtPerf=CxtPerf(~isnan(CxtPerf.session),:);

for c=1:size(CxtPerf,1)
CxtPerf.acc(c) = 2-nanmean(Behav.correctness(Behav.session==CxtPerf.session(c) & Behav.context==CxtPerf.context(c)));
end


%%
figure;

boxplot(CxtPerf.acc,CxtPerf.context)
hold on
scatter(CxtPerf.context,CxtPerf.acc)

anova(CxtPerf,'acc')

[h,p]=ttest2(CxtPerf.acc(CxtPerf.context==3),CxtPerf.acc(CxtPerf.context==1))