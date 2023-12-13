Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed];
ROOT.Rip = [ROOT.Save '\ripples_mat\R2'];
ROOT.Unit = [ROOT.Save '\units_mat\U0'];


% Rip_CA1 = readtable([ROOT.Save '\RipplesTable_' 'CA1_field' '_RDIs_UV_cell_HeteroIn_AllPopul.xlsx']);
% Rip_SUB = readtable([ROOT.Save '\RipplesTable_' 'SUB_field' '_RDIs_UV_cell_HeteroIn_AllPopul.xlsx']);
suff = '_forAnalysis'; 
% suff='_Behav';
% Rip_CA1 = readtable([ROOT.Save '\RipplesTable_' 'CA1' suff '.xlsx']);
% Rip_SUB = readtable([ROOT.Save '\RipplesTable_' 'SUB' suff '.xlsx']);
Rip_CA1 = readtable([ROOT.Rip '\RipplesTable_Behav_' 'CA1_speed_filtered' '.xlsx']);
Rip_SUB = readtable([ROOT.Rip '\RipplesTable_Behav_' 'SUB_speed_filtered' '.xlsx']);


pt = table;

Rip_CA1.SpikePerCell = Rip_CA1.spike./Rip_CA1.ensemble;
Rip_SUB.SpikePerCell = Rip_SUB.spike./Rip_SUB.ensemble;


Rip_CA1(Rip_CA1.rat==232 & Rip_CA1.session==4,:)=[];
Rip_SUB(Rip_SUB.rat==232 & Rip_SUB.session==4,:)=[];
% Rip_CA1 = Rip_CA1(Rip_CA1.NormAmp_v3>0,:);


%%
vrList = {'RippleDuration','MeanRaw','MeanFreq','RipplePower','spike','ensemble','SpikePerCell'};
for v=1:numel(vrList)
vr  = vrList{v};

pt.(vr)(1) = nanmean(Rip_SUB.(vr));
pt.(vr)(2) = nanstd(Rip_SUB.(vr));


pt.(vr)(3) = nanmean(Rip_CA1.(vr));
pt.(vr)(4) = nanstd(Rip_CA1.(vr));

x1 = Rip_SUB.(vr); x1(isnan(x1))=[];
x2 = Rip_CA1.(vr); x2(isnan(x2))=[];
[h,p] = ttest2(x1,x2);

pt.(vr)(5) = p;
end
% sum(Rip_CA1.nFields>=5)

%%
pt.RipplePower = pt.RipplePower/(10^-12)

%%

%% Ripple Property Bar graphs
for v=1:1
vr  = vrList{v};

x = [1 2];
data = [pt.(vr)(1); pt.(vr)(3)];
err = [pt.(vr)(2)/sqrt(size(Rip_SUB,1)); pt.(vr)(4)/sqrt(size(Rip_CA1,1))];

% b = bar(x,data); 
% errorbar(data,err,"LineStyle","none",'color','k')
% 
% b.FaceColor = "flat";
% b.CData(1,:) = CList(1,:);
% b.CData(2,:) = CList(2,:);
figure; hold on
% b=violinplot([ ;Rip_CA1.(vr)],[ones(size(Rip_SUB,1),1);2*ones(size(Rip_CA1,1),1)])
h2 = histogram(Rip_CA1.(vr),'binwidth',0.01,'Normalization','probability');
h = histogram(Rip_SUB.(vr),'binwidth',h2.BinWidth, 'binlimits',h2.BinLimits,'Normalization','probability');
legend({'CA1','SUB'})
 xlabel(vr)

[h,k] = kstest2(Rip_SUB.(vr),Rip_CA1.(vr));
[p,h,stats] = ranksum(Rip_SUB.(vr),Rip_CA1.(vr));
title ([vr '-Z=' jjnum2str(stats.zval,2) ',w=' jjnum2str(p,5) ',k=' jjnum2str(k,4)])

ylabel('proportion')

% ylim([180 200])
if ~exist([ROOT.Save '\Manuscript\Table1']), mkdir([ROOT.Save '\Manuscript\Table1']); end
xlim([0 0.2])
saveas(gca,[ROOT.Save '\Manuscript\Table1\(hist)(filtered)' vr '.png'])
saveas(gca,[ROOT.Save '\Manuscript\Table1\(hist)(filtered)' vr '.svg'])


vr
end
%% Ripple Rate_rat
RatList = unique([Rip_SUB.rat;Rip_CA1.rat]);
Rip_Rats=[];

for rid = 1:size(RatList,1)
    Rip_Rats(rid,1) = sum(Rip_SUB.rat==RatList(rid));
     Rip_Rats(rid,2) = sum(Rip_CA1.rat==RatList(rid));
end

boxplot(Rip_Rats);
hold on
scatter([1 1 1 1 1 2 2 2 2 2]',[Rip_Rats(:,1);Rip_Rats(:,2)])

x=Rip_Rats(:,1); y=Rip_Rats(:,2)

[p,h] = ranksum(x,y)
%% Ripple Rate_session
RatList = unique([Rip_SUB.rat;Rip_CA1.rat]);
SessionList = unique([Rip_SUB.session;Rip_CA1.session]);
Rip_Rats=[];
r=1;
for rid = 1:size(RatList,1)
    for sid=1:size(SessionList,1)
    Rip_Rats(r,1) = sum(Rip_SUB.rat==RatList(rid) & Rip_SUB.session==SessionList(sid));
     Rip_Rats(r,2) = sum(Rip_CA1.rat==RatList(rid) & Rip_CA1.session==SessionList(sid));
     Rip_Rats(r,3) = RatList(rid);
     Rip_Rats(r,4) = SessionList(sid);
     r=r+1;
    end
end

% Rip_Rats(sum(Rip_Rats,2)==0,:)=[];
x=Rip_Rats(:,1); y=Rip_Rats(:,2);
x(x==0)=nan; y(y==0)=nan;
figure; hold on
boxplot([x,y]);
scatter(ones(length(x),1),x,40,Rip_Rats(:,3),'filled')
scatter(ones(length(x),1)*2,y,40,Rip_Rats(:,3),'filled')
% ylim([0 450])

[p,h,stats] = ranksum(x,y)
legend
title (['ripples' '-Z=' jjnum2str(stats.zval,2) ',w=' jjnum2str(p,4)])

