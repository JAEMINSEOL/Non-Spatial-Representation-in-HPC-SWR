Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed ''];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Save '\ripples_mat\R4'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R11_sub'];
ROOT.Units = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];


thisRegion = 'CA1';
thisRegion2 = [thisRegion '_field'];
% RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_RDIs_UV_cell_HeteroIn_AllPopul' '.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable.CA1_field  = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);

thisRegion = 'SUB';
thisRegion2 = [thisRegion '_field'];
% RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_RDIs_UV_cell_HeteroIn_AllPopul' '.xlsx']);
UnitsTable.SUB = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable.SUB_field = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);

CList = [[139 193 69]/255; [29 111 169]/255];
%%
thisRegion = 'SUB'; cl=1;

thisRegion2 = [thisRegion '_field'];
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);

RT = RipplesTable.(thisRegion);
for uid=1:size(UnitsTable.(thisRegion),1)
    
    thisRips = ReactTable(find(cellfun(Params.cellfind(UnitsTable.(thisRegion).ID{uid}),ReactTable.UnitID)),[1 2]);

    if ~isempty(thisRips)
        for rid=1:size(thisRips,1)
            thisRips.RippleDuration(rid) = RT.RippleDuration(strcmp(RT.ID,thisRips.RippleID{rid}));
        end
        thisR =unique(thisRips);
    UnitsTable.(thisRegion).RipAvgFR(uid) = size(thisRips,1)/sum(thisR.RippleDuration);
    else
    UnitsTable.(thisRegion).RipAvgFR(uid) = 0;
    end
end
%%
thisRegion = 'SUB'; 
UT = UnitsTable.(thisRegion); UT(UT.RipAvgFR==0,:)=[];
scatter(UT.AvgFR,UT.meanSpks)
mdl = fitlm(UT.onMazeAvgFR,UT.RipAvgFR)
mdl.Rsquared.Adjusted
%% 각 unit별 전체 spike 수 계산
thisRegion = 'CA1'; cl=1;
RT = RipplesTable.(thisRegion);
thisRegion2 = [thisRegion '_field'];
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);
ReactSpks=[];

for uid=1:size(UnitsTable.(thisRegion),1)
    thisReact = ReactTable(find(cellfun(Params.cellfind(UnitsTable.(thisRegion).ID{uid}),ReactTable.UnitID)),[1 2]);

    thisR= unique(thisReact.RippleID);
    if ~isempty(thisR)
        ids=[];
        for rid=1:size(thisR,1)
            if RT.pRDI_L_UV(strcmp(RT.ID,thisR{rid}))<0.05
                ids(rid)=1;
            else, ids(rid)=0;
            end
        end
        if sum(ids)>0
        thisR=thisR(find(ids),:);
    for rct = 1:size(thisR,1)
         n = size(thisReact(find(cellfun(Params.cellfind(thisR{rct}), thisReact.RippleID)),1),1);
         f = UnitsTable.(thisRegion).NumField(uid);
ReactSpks = [ReactSpks; [n f 1]];

    end
        end
    end
end

thisRegion = 'SUB'; cl=1;
RT = RipplesTable.(thisRegion);
thisRegion2 = [thisRegion '_field'];
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);

for uid=1:size(UnitsTable.(thisRegion),1)
    thisReact = ReactTable(find(cellfun(Params.cellfind(UnitsTable.(thisRegion).ID{uid}),ReactTable.UnitID)),[1 2]);

    thisR= unique(thisReact.RippleID);
    if ~isempty(thisR)
            ids=[];
        for rid=1:size(thisR,1)
            if RT.pRDI_L_UV(strcmp(RT.ID,thisR{rid}))<0.05
                ids(rid)=1;
            else, ids(rid)=0;
            end
        end
                if sum(ids)>0
        thisR=thisR(find(ids),:);
    for rct = 1:size(thisR,1)
         n = size(thisReact(find(cellfun(Params.cellfind(thisR{rct}), thisReact.RippleID)),1),1);
         f = UnitsTable.(thisRegion).NumField(uid);
ReactSpks = [ReactSpks; [n f 0]];

    end
                end
    end
end

%% 각 unit별 전체 spike 계산_plot
ReactSpks = ReactSpks(ReactSpks(:,2)>0,:);
RS = ReactSpks;

subplot(2,2,1); hold on
RSs = RS(RS(:,3)==0,:);
x1 = RSs(:,2); y1=RSs(:,1);
t1 = [mean(y1(x1==1)) mean(y1(x1==2)) mean(y1(x1>=3))];
e1 = [std(y1(x1==1))/sqrt(sum(y1(x1==1))) std(y1(x1==2))/sqrt(sum(y1(x1==2))) std(y1(x1>=3))/sqrt(sum(y1(x1>=3)))];
b = bar([1 2 3], t1);
b.FaceColor = 'flat';
b.CData = CList(1,:);
xpos = b(1).XData + b(1).XOffset;
errorbar(xpos, t1, e1, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
set(gca,'fontsize',12,'fontweight','b')
xticks([1 2 3]);
xticklabels({'1' '2' '>3'})
title('SUB'); xlabel('# fields'); ylim([0 2]);ylabel('spikes in a ripple')

subplot(2,2,3); hold on
RSs = RS(RS(:,3)==0,:);
x1 = RSs(:,2); y1=RSs(:,1);
t1 = [mean(y1(x1==1)) mean(y1(x1>=2)) ];
e1 = [std(y1(x1==1))/sqrt(sum(y1(x1==1))) std(y1(x1>=2))/sqrt(sum(y1(x1>=2)))];
b = bar([1 2], t1);
b.FaceColor = 'flat';
b.CData = CList(1,:);
xpos = b(1).XData + b(1).XOffset;
errorbar(xpos, t1, e1, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
set(gca,'fontsize',12,'fontweight','b')
xticks([1 2]);
xticklabels({'SF' 'MF'})
title('SUB'); xlabel('# fields'); ylim([0 2]);ylabel('spikes in a ripple')

subplot(2,2,2); hold on
RSc = RS(RS(:,3)==1,:);
x2 = RSc(:,2); y2=RSc(:,1);
t1 = [mean(y2(x2==1)) mean(y2(x2==2)) mean(y2(x2>=3))];
e1 = [std(y2(x2==1))/sqrt(sum(y2(x2==1))) std(y2(x2==2))/sqrt(sum(y2(x2==2))) std(y2(x2>=3))/sqrt(sum(y2(x2>=3)))];
b = bar([1 2 3], t1);
b.FaceColor = 'flat';
b.CData = CList(2,:);
xpos = b(1).XData + b(1).XOffset;
errorbar(xpos, t1, e1, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
set(gca,'fontsize',12,'fontweight','b')
xticks([1 2 3]);
xticklabels({'1' '2' '>3'})
title('CA1'); xlabel('Unit type'); ylim([0 2]);ylabel('spikes in a ripple')

subplot(2,2,4); hold on
RSc = RS(RS(:,3)==1,:);
x2 = RSc(:,2); y2=RSc(:,1);
t1 = [mean(y2(x2==1)) mean(y2(x2>=2)) ];
e1 = [std(y2(x2==1))/sqrt(sum(y2(x2==1))) std(y2(x2>=2))/sqrt(sum(y2(x2>=2)))];
b = bar([1 2], t1);
b.FaceColor = 'flat';
b.CData = CList(2,:);
xpos = b(1).XData + b(1).XOffset;
errorbar(xpos, t1, e1, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
set(gca,'fontsize',12,'fontweight','b')
xticks([1 2]);
xticklabels({'SF' 'MF'})
title('CA1'); xlabel('Unit type'); ylim([0 2]);ylabel('spikes in a ripple')

sgtitle('NSp React (Left)')
%% 각 unit별 평균 spike 수 계산
thisRegion = 'CA1'; cl=1;
RT = RipplesTable.(thisRegion);
thisRegion2 = [thisRegion '_field'];
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);

RT = RipplesTable.(thisRegion);
for uid=1:size(UnitsTable.(thisRegion),1)
    
    thisRips = ReactTable(find(cellfun(Params.cellfind(UnitsTable.(thisRegion).ID{uid}),ReactTable.UnitID)),[1 2]);

    if ~isempty(thisRips)
    UnitsTable.(thisRegion).meanSpks(uid) = size(thisRips,1)/size(unique(thisRips.RippleID),1);
    else
    UnitsTable.(thisRegion).meanSpks(uid) = 0;
    end
end



%% 각 unit별 평균 spike 수 계산_plot

UT = UnitsTable.SUB; 
UT(~UT.meanSpks | ~UT.NumField,:)=[];
x1 = UT.NumField; y1 = UT.meanSpks;

figure('Position',[2400,-500,560,1200])
subplot(3,2,1)
scatter(x1,y1,20,'k','filled')
title('SUB'); xlim([0 6]); ylim([0 4]); xticks([1:5]); xlabel('# fields'); ylabel('mean spikes in ripples')
set(gca,'fontsize',12,'fontweight','b')

UT = UnitsTable.CA1; 
UT(~UT.meanSpks | ~UT.NumField,:)=[];
x2 = UT.NumField; y2 = UT.meanSpks;

subplot(3,2,2)
scatter(x2,y2,20,'k','filled')
title('CA1'); xlim([0 6]); ylim([0 4]);xticks([1:5]); xlabel('# fields'); ylabel('mean spikes in ripples')
set(gca,'fontsize',12,'fontweight','b')

subplot(3,2,3); hold on
t1 = [mean(y1(x1==1)) mean(y1(x1==2)) mean(y1(x1>=3))];
e1 = [std(y1(x1==1))/sqrt(sum(y1(x1==1))) std(y1(x1==2))/sqrt(sum(y1(x1==2))) std(y1(x1>=3))/sqrt(sum(y1(x1>=3)))];
b = bar([1 2 3], t1);
b.FaceColor = 'flat';
b.CData = CList(1,:);
xpos = b(1).XData + b(1).XOffset;
errorbar(xpos, t1, e1, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
set(gca,'fontsize',12,'fontweight','b')
xticks([1 2 3]);
xticklabels({'1' '2' '>3'})
title('SUB'); xlabel('# fields'); ylim([0 3]);ylabel('mean spikes in ripples')
[h,p] = ranksum(y2(x2>=3),y2(x2==1))

subplot(3,2,4); hold on
t1 = [mean(y2(x2==1)) mean(y2(x2==2)) mean(y2(x2>=3))];
e1 = [std(y2(x1==1))/sqrt(sum(y1(x1==1))) std(y1(x1==2))/sqrt(sum(y1(x1==2))) std(y1(x1>=3))/sqrt(sum(y1(x1>=3)))];
b = bar([1 2 3], t1);
b.FaceColor = 'flat';
b.CData = CList(2,:);
xpos = b(1).XData + b(1).XOffset;
errorbar(xpos, t1, e1, 'LineStyle', 'none', 'Color', 'k', 'LineWidth', 1);
set(gca,'fontsize',12,'fontweight','b')
xticks([1 2 3]);
xticklabels({'1' '2' '>3'})
title('CA1'); xlabel('# fields');  ylim([0 3]);ylabel('mean spikes in ripples')

subplot(3,2,5)
boxplot(y1,x1>1)
xticks([1 2]); xlim([.5 2.5])
xticklabels({'SF' 'MF'})
title('SUB'); xlabel('cell type');  ylim([0 4]);ylabel('mean spikes in ripples')
set(gca,'fontsize',12,'fontweight','b')
[p,h] = ranksum(y1(x1==1), y1(x1>1));
text(1.2,3,['p= ' jjnum2str(p,4)])

subplot(3,2,6); hold on
boxplot(y2,x2>2)
xticks([1 2]); xlim([.5 2.5])
xticklabels({'SF' 'MF'})
title('CA1'); xlabel('cell type');  ylim([0 4]);ylabel('mean spikes in ripples')
set(gca,'fontsize',12,'fontweight','b')
[p,h] = ranksum(y2(x2==1), y2(x2>1));
text(1.2,3,['p= ' jjnum2str(p,4)])


%% 각 unit pair의 co-firing 분석 
thisRegion = 'SUB'; cl=0;
RT = RipplesTable.(thisRegion);
thisRegion2 = [thisRegion '_field'];
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion2 '.xlsx']);
UT = UnitsTable.(thisRegion2); 
RT = RipplesTable.(thisRegion);

UnitPair=table; u=1;
for uid=1:size(UT,1)
    p0 = sum(RT.rat==UT.rat(uid) & RT.session==UT.session(uid));
    thisUnits = UT(UT.rat==UT.rat(uid) & UT.session==UT.session(uid),:);
    thisRips = unique(ReactTable(find(cellfun(Params.cellfind(UT.ID{uid}),ReactTable.UnitID)),[1 2]));
    p1 = size(thisRips,1);
    for uid2=1:size(thisUnits,1)
        if strncmp(thisUnits.ID{uid2},UT.ID{uid},12), continue; end

        thisRips2=unique(ReactTable(find(cellfun(Params.cellfind(thisUnits.ID{uid2}),ReactTable.UnitID)),[1 2]));
        CoActRips = intersect(thisRips.RippleID,thisRips2.RippleID);

        q1 = size(thisRips2,1);
        pq = size(CoActRips,1);

        UnitPair.UID1{u}=UT.ID{uid};
        UnitPair.UID2{u}=thisUnits.ID{uid2};
        UnitPair.p(u)=p1;
        UnitPair.q(u)=q1;
        UnitPair.pq(u)=pq;
        UnitPair.un(u) = size(union(thisRips.RippleID,thisRips2.RippleID),1);
        UnitPair.p0(u)=p0;

        UnitPair.L1(u)=UT.RDI_LScene(uid);
        UnitPair.R1(u)=UT.RDI_RScene(uid);
        UnitPair.C1(u)=UT.RDI_LR(uid);

          UnitPair.L2(u)=thisUnits.RDI_LScene(uid2);
        UnitPair.R2(u)=thisUnits.RDI_RScene(uid2);
        UnitPair.C2(u)=thisUnits.RDI_LR(uid2);

        UnitPair.region(u) = cl;
             UnitPair.NumField1(u) = UT.NumField(uid);
             UnitPair.NumField2(u) = thisUnits.NumField(uid2);

%         UnitPair.NumField1(u) = sum(UT.AvgFR(uid)==UT.AvgFR);
%          UnitPair.NumField2(u) = sum(thisUnits.AvgFR(uid2)==thisUnits.AvgFR);
        
        u=u+1;
    end
end
writetable(UnitPair,[ROOT.Save '\UnitPair_' thisRegion2 '.xlsx'],'writemode','replacefile');
%%
thisRegion = 'CA1'; cl=0;
RT = RipplesTable.(thisRegion);
thisRegion2 = [thisRegion '_field'];
UnitPair_ca1 =  readtable([ROOT.Save '\UnitPair_' thisRegion2 '.xlsx']);
UnitPair_ca12 =  readtable([ROOT.Save '\UnitPair_' thisRegion '.xlsx']);

thisRegion = 'SUB'; cl=0;
RT = RipplesTable.(thisRegion);
thisRegion2 = [thisRegion '_field'];
UnitPair_sub =  readtable([ROOT.Save '\UnitPair_' thisRegion2 '.xlsx']);
UnitPair_sub2 =  readtable([ROOT.Save '\UnitPair_' thisRegion '.xlsx']);

%% Plot_Unit Co-firing
UnitPair=[UnitPair_sub];
upl = UnitPair;
upl = upl(min([abs(upl.L1),abs(upl.L2)],[],2)>=0.1,:);

d1l = [upl((upl.L1 .* upl.L2) <0,:)]; 
d2l = [upl((upl.L1 .* upl.L2) >0,:)]; 

upl = UnitPair;
upl = upl(min([abs(upl.R1),abs(upl.R2)],[],2)>=0.1,:);
d1r = [upl((upl.R1 .* upl.R2) <0,:)]; 
d2r = [upl((upl.R1 .* upl.R2) >0,:)]; 

upl = UnitPair;
upl = upl(min([abs(upl.C1),abs(upl.C2)],[],2)>=0.1,:);
d1c = [upl((upl.C1 .* upl.C2) <0,:)]; 
d2c = [upl((upl.C1 .* upl.C2) >0,:)]; 

d1 = [d1l; d1r; d1c];
d2 = [d2l; d2r; d2c];

    x1 = d1.pq./min([d1.p,d1.q],[],2);
x2 = d2.pq./min([d2.p,d2.q],[],2);

x1(x1==1|x1==0)=[];
x2(x2==1|x2==0)=[];

% y = abs(upl.L1-upl.L2);
% x=upl.pq./upl.un;
% id=(isnan(x)|isnan(y)|x==1|x==0);
% x(id)=[]; y(id)=[];
% scatter(y,x)
% corrcoef(x,y)
figure; 
subplot(1,2,1); hold on
data = [nanmean(x1) nanmean(x2)];
err = [nanstd(x1)/sqrt(size(x1,1)) nanstd(x2)/sqrt(size(x2,1)) ];
b = bar(data);

xpos = b.XData + b.XOffset;


    errorbar(xpos, data, err, 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);


[h,p] = ttest2(x1,x2)
xticks([1 2]); xticklabels({'Hetero-pair','Homo-pair'})
ylabel('Co-firing probability')
set(gca,'fontsize',12,'fontweight','b')
ylim([0 .5])
title('SUB')



UnitPair=[UnitPair_ca1];
upl = UnitPair;
upl = upl(min([abs(upl.L1),abs(upl.L2)],[],2)>=0.1,:);

d1l = [upl((upl.L1 .* upl.L2) <0,:)]; 
d2l = [upl((upl.L1 .* upl.L2) >0,:)]; 

upl = UnitPair;
upl = upl(min([abs(upl.R1),abs(upl.R2)],[],2)>=0.1,:);
d1r = [upl((upl.R1 .* upl.R2) <0,:)]; 
d2r = [upl((upl.R1 .* upl.R2) >0,:)]; 

upl = UnitPair;
upl = upl(min([abs(upl.C1),abs(upl.C2)],[],2)>=0.1,:);
d1c = [upl((upl.C1 .* upl.C2) <0,:)]; 
d2c = [upl((upl.C1 .* upl.C2) >0,:)]; 

d1 = [d1l; d1r; d1c];
d2 = [d2l; d2r; d2c];

   x1 = d1.pq./min([d1.p,d1.q],[],2);
x2 = d2.pq./min([d2.p,d2.q],[],2);

x1(x1==1|x1==0)=[];
x2(x2==1|x2==0)=[];

% y = abs(upl.L1-upl.L2);
% x=upl.pq./upl.un;
% id=(isnan(x)|isnan(y)|x==1|x==0);
% x(id)=[]; y(id)=[];
% scatter(y,x)
% corrcoef(x,y)
subplot(1,2,2); hold on
data = [nanmean(x1) nanmean(x2)];
err = [nanstd(x1)/sqrt(size(x1,1)) nanstd(x2)/sqrt(size(x2,1)) ];
b = bar(data);

xpos = b.XData + b.XOffset;


    errorbar(xpos, data, err, 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);


[h,p] = ttest2(x1,x2)
xticks([1 2]); xticklabels({'Hetero-pair','Homo-pair'})
ylabel('Co-firing probability')
set(gca,'fontsize',12,'fontweight','b')
ylim([0 .5])
title('CA1')
%% Plot_Unit Co-firing_MFSF
UnitPair=[UnitPair_sub2];
upl = UnitPair;
upl = upl(min([abs(upl.L1),abs(upl.L2)],[],2)>=0.1,:);

d1l = [upl(upl.NumField1==1 & upl.NumField2==1,:)]; 
d2l = [upl(upl.NumField1>1 & upl.NumField2>1,:)]; 
d3 = [upl((upl.NumField1>1 & upl.NumField2==1) | upl.NumField2>1 & upl.NumField1==1,:)]; 

d1 = [d1l];
d2 = [d2l];

    x1 = d1.pq./min([d1.p,d1.q],[],2);
x2 = d2.pq./min([d2.p,d2.q],[],2);
x3 = d3.pq./min([d3.p,d3.q],[],2);

x1(x1==1|x1==0)=[];
x2(x2==1|x2==0)=[];
x3(x3==1|x3==0)=[];


figure; 
subplot(1,2,1); hold on
data = [nanmean(x1) nanmean(x2) nanmean(x3)];
err = [nanstd(x1)/sqrt(size(x1,1)) nanstd(x2)/sqrt(size(x2,1)) nanstd(x3)/sqrt(size(x3,1))];
b = bar(data);

xpos = b.XData + b.XOffset;


    errorbar(xpos, data, err, 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);


[h,p] = ttest2(x1,x2)
xticks([1 2 3]); xticklabels({'SF-SF','MF-MF','SF-MF'})
ylabel('Co-firing probability')
set(gca,'fontsize',12,'fontweight','b')
ylim([0 .6])
title('SUB')

UnitPair=[UnitPair_ca12];
upl = UnitPair;

d1l = [upl(upl.NumField1==1 & upl.NumField2==1,:)]; 
d2l = [upl(upl.NumField1>1 & upl.NumField2>1,:)]; 
d3 = [upl((upl.NumField1>1 & upl.NumField2==1) | upl.NumField2>1 & upl.NumField1==1,:)]; 

d1 = [d1l];
d2 = [d2l];

    x1 = d1.pq./min([d1.p,d1.q],[],2);
x2 = d2.pq./min([d2.p,d2.q],[],2);
x3 = d3.pq./min([d3.p,d3.q],[],2);

x1(x1==1|x1==0)=[];
x2(x2==1|x2==0)=[];
x3(x3==1|x3==0)=[];


subplot(1,2,2); hold on
data = [nanmean(x1) nanmean(x2) nanmean(x3)];
err = [nanstd(x1)/sqrt(size(x1,1)) nanstd(x2)/sqrt(size(x2,1)) nanstd(x3)/sqrt(size(x3,1))];
b = bar(data);

xpos = b.XData + b.XOffset;


    errorbar(xpos, data, err, 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);


[h,p] = ttest2(x1,x2)
xticks([1 2 3]); xticklabels({'SF-SF','MF-MF','SF-MF'})
ylabel('Co-firing probability')
set(gca,'fontsize',12,'fontweight','b')
ylim([0 .6])
title('CA1')


%%
t1 = UnitsTable_sub_f; t2 = UnitsTable_ca1_f;

[p1s,h]=ranksum(t1.RDI_LScene(abs(t1.RDI_LScene)>.1),0);
[p2s,h]=ranksum(t1.RDI_RScene(abs(t1.RDI_RScene)>.1),0);
[p3s,h]=ranksum(t1.RDI_LR(abs(t1.RDI_LR)>.1),0);

[p1c,h]=ranksum(t2.RDI_LScene(abs(t2.RDI_LScene)>.1),0);
[p2c,h]=ranksum(t2.RDI_RScene(abs(t2.RDI_RScene)>.1),0);
[p3c,h]=ranksum(t2.RDI_LR(abs(t2.RDI_LR)>.1),0);

figure('position',[2088,206,1385,612])
subplot(2,3,1)
histogram(t1.RDI_LScene(abs(t1.RDI_LScene)>.1),'binwidth',0.1,'facecolor',CList(1,:),'Normalization','probability')
xlim([-1.2 1.2]); xlabel('Left scene selectivity'); ylim([0 .2])
title('unit selectivity distribution, Subiculum')
subplot(2,3,2)
histogram(t1.RDI_RScene(abs(t1.RDI_RScene)>.1),'binwidth',0.1,'facecolor',CList(1,:),'Normalization','probability')
xlim([-1.2 1.2]); xlabel('Right scene selectivity'); ylim([0 .2])
title('unit selectivity distribution, Subiculum')
subplot(2,3,3)
histogram(t1.RDI_LR(abs(t1.RDI_LR)>.1),'binwidth',0.1,'facecolor',CList(1,:),'Normalization','probability')
xlim([-1.2 1.2]); xlabel('Choice selectivity'); ylim([0 .2])
title('unit selectivity distribution, Subiculum')
subplot(2,3,4)
histogram(t2.RDI_LScene(abs(t2.RDI_LScene)>.1),'binwidth',0.1,'facecolor',CList(2,:),'Normalization','probability')
xlim([-1.2 1.2]); xlabel('Left scene selectivity'); ylim([0 .2])
title('unit selectivity distribution, CA1')
subplot(2,3,5)
histogram(t2.RDI_RScene(abs(t2.RDI_RScene)>.1),'binwidth',0.1,'facecolor',CList(2,:),'Normalization','probability')
xlim([-1.2 1.2]); xlabel('Right scene selectivity'); ylim([0 .2])
title('unit selectivity distribution, CA1')
subplot(2,3,6)
histogram(t2.RDI_LR(abs(t2.RDI_LR)>.1),'binwidth',0.1,'facecolor',CList(2,:),'Normalization','probability')
xlim([-1.2 1.2]); xlabel('Choice selectivity'); ylim([0 .2])
title('unit selectivity distribution, CA1')

%%
for uid = 1:size(UnitsTable_ca1)
thisUnit = UnitsTable_ca1(uid,:);

UnitsTable_ca1.NumField(uid) =...
    size(find(thisUnit.nSpks==UnitsTable_ca1_f.nSpks & thisUnit.SpkWidth==UnitsTable_ca1_f.SpkWidth),1);
end


UnitsTable_ca1a = UnitsTable_ca1(UnitsTable_ca1.NumField>0,:);