Initial_SWRFilter_common;
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
Behavior = readtable([ROOT.Behav '/BehavTable.xlsx']);
%%
for t = 1:size(Behavior)
    if ismember(str2double(Behavior.TrialID{t}(1:3)),[232 234 295 561])
        Behavior.rat(t) = str2double(Behavior.TrialID{t}(1:3));
        Behavior.session(t) = str2double(Behavior.TrialID{t}(5:6));
    Behavior.RS(t) = str2double(Behavior.TrialID{t}(1:3))*1e2+str2double(Behavior.TrialID{t}(5:6));
    else
        Behavior.RS(t) =nan;
    end
end

%%

CxtPerf = table;
CxtPerf = unique(Behavior(:,[17,2]),'rows');

CxtPerf=CxtPerf(~isnan(CxtPerf.RS),:);


for c=1:size(CxtPerf,1)
    CxtPerf.rat(c) = Behavior.rat(find(Behavior.RS==CxtPerf.RS(c),1));
    CxtPerf.session(c) = Behavior.session(find(Behavior.RS==CxtPerf.RS(c),1));
CxtPerf.acc(c) = 2-nanmean(Behavior.correctness(Behavior.RS==CxtPerf.RS(c) & Behavior.context==CxtPerf.context(c)));

CxtPerf.mduration_ITI(c) = nanmean(Behavior.duration_ITI(Behavior.RS==CxtPerf.RS(c) & Behavior.context==CxtPerf.context(c) & Behavior.correctness==1));
CxtPerf.sduration_ITI(c) = nanstd(Behavior.duration_ITI(Behavior.RS==CxtPerf.RS(c) & Behavior.context==CxtPerf.context(c) & Behavior.correctness==1));
end

for t = 1:size(Behavior,1)
    if ismember(str2double(Behavior.TrialID{t}(1:3)),[232 234 295 561])
        c = find(Behavior.RS(t)==CxtPerf.RS);
        if Behavior.duration_ITI(t)>mean(CxtPerf.mduration_ITI(c))+mean(CxtPerf.sduration_ITI(c))
            Behavior.duration_ITI(t)=nan;
        end
    else

    end
end

for c=1:size(CxtPerf,1)
CxtPerf.mduration_ITI(c) = nanmean(Behavior.duration_ITI(Behavior.RS==CxtPerf.RS(c) & Behavior.context==CxtPerf.context(c) & Behavior.correctness==1));
CxtPerf.sduration_ITI(c) = nanstd(Behavior.duration_ITI(Behavior.RS==CxtPerf.RS(c) & Behavior.context==CxtPerf.context(c) & Behavior.correctness==1));
end

CxtPerf.rat(CxtPerf.rat==234)=400;
%%
figure;

boxplot(CxtPerf.acc,CxtPerf.context)
hold on
scatter(CxtPerf.context,CxtPerf.acc,20,CxtPerf.rat,'filled')
%%
figure;

boxplot(CxtPerf.mduration_ITI,CxtPerf.context)
hold on
scatter(CxtPerf.context,CxtPerf.mduration_ITI,20,CxtPerf.rat,'filled')
ylim([5 20])


%%
id = Behavior.rat==232 & Behavior.session==4;
x = Behavior(id,[8,2]);
x(isnan(x.duration_ITI),:)=[];
anova(x,'duration_ITI')

% colormap([0 0 0; 1 0 0; 0 0 1; 0 1 0; 1 1 0;0 1 1])
