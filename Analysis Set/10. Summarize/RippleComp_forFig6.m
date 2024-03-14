%%
% RipplesTable_replay.SUB = readtable([ROOT.Save '\RipplesTable_' 'SUB_refCA1' '_field_RDIs_UV_pc.xlsx']);
% RipplesTable_replay.CA1 = readtable([ROOT.Save '\RipplesTable_' 'CA1' '_field_RDIs_UV_pc.xlsx']);

addpath('D:\HPC-SWR project\Analysis Program')
Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed];

%%
clear RipplesTable
RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_' 'SUB' '_forAnalysis.xlsx']);
RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_' 'CA1' '_forAnalysis.xlsx']);


Rip_SUB = RipplesTable.SUB;
Rip_CA1 = RipplesTable.CA1;

Rip_SUB(Rip_SUB.rat==232 & Rip_SUB.session==4,:) = [];
Rip_CA1(Rip_CA1.rat==232 & Rip_CA1.session==4,:) = [];
al = 0.05;
%%

Rip_SUB.pRDI_L_max(Rip_SUB.nRDI_L_max<5) = nan;
Rip_SUB.pRDI_R_max(Rip_SUB.nRDI_R_max<5) = nan;
Rip_SUB.pRDI_C_max(Rip_SUB.nRDI_C_max<5) = nan;

Rip_CA1.pRDI_L_max(Rip_CA1.nRDI_L_max<5) = nan;
Rip_CA1.pRDI_R_max(Rip_CA1.nRDI_R_max<5) = nan;
Rip_CA1.pRDI_C_max(Rip_CA1.nRDI_C_max<5) = nan;
%%
Rip_SUB.pBinomDev_L_UV(Rip_SUB.nRDI_L_max<5) = nan;
Rip_SUB.pBinomDev_R_UV(Rip_SUB.nRDI_R_max<5) = nan;
Rip_SUB.pBinomDev_C_UV(Rip_SUB.nRDI_C_max<5) = nan;

Rip_CA1.pBinomDev_L_UV(Rip_CA1.nRDI_L_max<5) = nan;
Rip_CA1.pBinomDev_R_UV(Rip_CA1.nRDI_R_max<5) = nan;
Rip_CA1.pBinomDev_C_UV(Rip_CA1.nRDI_C_max<5) = nan;


%%
Rip_SUB.pRatio_L_UV(Rip_SUB.nRDI_L_max<5) = nan;
Rip_SUB.pRatio_R_UV(Rip_SUB.nRDI_R_max<5) = nan;
Rip_SUB.pRatio_C_UV(Rip_SUB.nRDI_C_max<5) = nan;

Rip_CA1.pRatio_L_UV(Rip_CA1.nRDI_L_max<5) = nan;
Rip_CA1.pRatio_R_UV(Rip_CA1.nRDI_R_max<5) = nan;
Rip_CA1.pRatio_C_UV(Rip_CA1.nRDI_C_max<5) = nan;
% writetable(RTsub,[ROOT.Save '\RipplesTable_SUB_refCA1_forAnalysis' '.xlsx'],'Writemode','replacefile');
% writetable(RTca1,[ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx'],'Writemode','replacefile');
%%
rats = unique([Rip_SUB.rat;Rip_CA1.rat]);


Sel_dist_sub = table;

for rid = 1:size(rats,1)
    Sel_dist_sub.Rat(rid) = rats(rid);
    Sel_dist_sub.nRipples(rid) = sum(Rip_SUB.rat==rats(rid) );
%         Sel_dist_sub.nValid(rid) = sum(RTsub.rat==rats(rid) & (RTsub.nRDI_L_UV>=5 | RTsub.nRDI_R_UV>=5 | RTsub.nRDI_C_UV>=5));
Sel_dist_sub.LeftOnly(rid) = sum(~(Rip_SUB.DecodingP_all<al) & Rip_SUB.rat==rats(rid) & Rip_SUB.pRDI_L_UV<al & ~(Rip_SUB.pRDI_R_UV<al) & ~(Rip_SUB.pRDI_C_UV<al));
Sel_dist_sub.RightOnly(rid) = sum(~(Rip_SUB.DecodingP_all<al) & Rip_SUB.rat==rats(rid) & ~(Rip_SUB.pRDI_L_UV<al) & Rip_SUB.pRDI_R_UV<al & ~(Rip_SUB.pRDI_C_UV<al));
Sel_dist_sub.ChoiceOnly(rid) = sum(~(Rip_SUB.DecodingP_all<al) & Rip_SUB.rat==rats(rid) & ~(Rip_SUB.pRDI_L_UV<al) & ~(Rip_SUB.pRDI_R_UV<al) & Rip_SUB.pRDI_C_UV<al);

Sel_dist_sub.LeftChoice(rid) = sum(~(Rip_SUB.DecodingP_all<al) & Rip_SUB.rat==rats(rid) & Rip_SUB.pRDI_L_UV<al & ~(Rip_SUB.pRDI_R_UV<al) & Rip_SUB.pRDI_C_UV<al);
Sel_dist_sub.RightChoice(rid) = sum(~(Rip_SUB.DecodingP_all<al) & Rip_SUB.rat==rats(rid) & ~(Rip_SUB.pRDI_L_UV<al) & Rip_SUB.pRDI_R_UV<al & Rip_SUB.pRDI_C_UV<al);
Sel_dist_sub.LeftRight(rid) = sum(~(Rip_SUB.DecodingP_all<al) & Rip_SUB.rat==rats(rid) & Rip_SUB.pRDI_L_UV<al & Rip_SUB.pRDI_R_UV<al & ~(Rip_SUB.pRDI_C_UV<al));

Sel_dist_sub.LRC(rid) = sum(~(Rip_SUB.DecodingP_all<al) & Rip_SUB.rat==rats(rid) & Rip_SUB.pRDI_L_UV<al & Rip_SUB.pRDI_R_UV<al & Rip_SUB.pRDI_C_UV<al);

Sel_dist_sub.Left_Z(rid) = sum(~(Rip_SUB.DecodingP_all<al) & Rip_SUB.mRDI_L_UV>0 & Rip_SUB.rat==rats(rid) & Rip_SUB.pRDI_L_UV<al & ~(Rip_SUB.pRDI_R_UV<al) & ~(Rip_SUB.pRDI_C_UV<al));
Sel_dist_sub.Left_B(rid) = sum(~(Rip_SUB.DecodingP_all<al) & Rip_SUB.mRDI_L_UV<0 & Rip_SUB.rat==rats(rid) & Rip_SUB.pRDI_L_UV<al & ~(Rip_SUB.pRDI_R_UV<al) & ~(Rip_SUB.pRDI_C_UV<al));


end





Sel_dist_ca1 = table;

for rid = 1:size(rats,1)
    Sel_dist_ca1.Rat(rid) = rats(rid);
    Sel_dist_ca1.nRipples(rid) = sum(Rip_CA1.rat==rats(rid));
%     Sel_dist_ca1.nValid(rid) = sum(RTca1.rat==rats(rid) & (RTca1.nRDI_L_UV>=5 | RTca1.nRDI_R_UV>=5 | RTca1.nRDI_C_UV>=5));
Sel_dist_ca1.LeftOnly(rid) = sum(~(Rip_CA1.DecodingP_all<al) & Rip_CA1.rat==rats(rid) & Rip_CA1.pRDI_L_UV<al & ~(Rip_CA1.pRDI_R_UV<al) & ~(Rip_CA1.pRDI_C_UV<al));
Sel_dist_ca1.RightOnly(rid) = sum(~(Rip_CA1.DecodingP_all<al) & Rip_CA1.rat==rats(rid) & ~(Rip_CA1.pRDI_L_UV<al) & Rip_CA1.pRDI_R_UV<al & ~(Rip_CA1.pRDI_C_UV<al));
Sel_dist_ca1.ChoiceOnly(rid) = sum(~(Rip_CA1.DecodingP_all<al) & Rip_CA1.rat==rats(rid) & ~(Rip_CA1.pRDI_L_UV<al) & ~(Rip_CA1.pRDI_R_UV<al) & Rip_CA1.pRDI_C_UV<al);

Sel_dist_ca1.LeftChoice(rid) = sum(~(Rip_CA1.DecodingP_all<al) & Rip_CA1.rat==rats(rid) & Rip_CA1.pRDI_L_UV<al & ~(Rip_CA1.pRDI_R_UV<al) & Rip_CA1.pRDI_C_UV<al);
Sel_dist_ca1.RightChoice(rid) = sum(~(Rip_CA1.DecodingP_all<al) & Rip_CA1.rat==rats(rid) & ~(Rip_CA1.pRDI_L_UV<al) & Rip_CA1.pRDI_R_UV<al & Rip_CA1.pRDI_C_UV<al);
Sel_dist_ca1.LeftRight(rid) = sum(~(Rip_CA1.DecodingP_all<al) & Rip_CA1.rat==rats(rid) & Rip_CA1.pRDI_L_UV<al & Rip_CA1.pRDI_R_UV<al & ~(Rip_CA1.pRDI_C_UV<al));

Sel_dist_ca1.LRC(rid) = sum(~(Rip_CA1.DecodingP_all<al) & Rip_CA1.rat==rats(rid) & Rip_CA1.pRDI_L_UV<al & Rip_CA1.pRDI_R_UV<al & Rip_CA1.pRDI_C_UV<al);
end
%%
T0 = Rip_SUB;
T1 = Rip_CA1;

%%
T0.mRDI_L_UV(T0.nFields<3)=0;
T1.mRDI_L_UV(T1.nFields<3)=0;
T0.mRDI_R_UV(T0.nFields<3)=0;
T1.mRDI_R_UV(T1.nFields<3)=0;
T0.mRDI_C_UV(T0.nFields<3)=0;
T1.mRDI_C_UV(T1.nFields<3)=0;

%%
T0.mRDI_L_UV(isnan(T0.mRDI_L_UV))=0;
T1.mRDI_L_UV(isnan(T1.mRDI_L_UV))=0;
T0.mRDI_R_UV(isnan(T0.mRDI_R_UV))=0;
T1.mRDI_R_UV(isnan(T1.mRDI_R_UV))=0;
T0.mRDI_C_UV(isnan(T0.mRDI_C_UV))=0;
T1.mRDI_C_UV(isnan(T1.mRDI_C_UV))=0;

%%
T0.mRDI_L_UV(isnan(T0.pRDI_L_UV))=nan;
T1.mRDI_L_UV(isnan(T1.pRDI_L_UV))=nan;
T0.mRDI_R_UV(isnan(T0.pRDI_R_UV))=nan;
T1.mRDI_R_UV(isnan(T1.pRDI_R_UV))=nan;
T0.mRDI_C_UV(isnan(T0.pRDI_C_UV))=nan;
T1.mRDI_C_UV(isnan(T1.pRDI_C_UV))=nan;

%%
crit='pBinom';
varList  = {'_L','_R', '_C'};
for v=1:3

var = varList{v};
% R0 = T0(T0.(['p' crit var])<al,:);
% R1 = T1(T1.(['p' crit var])<al,:);

T0.([crit var '_' crit2])(abs(T0.(['n_fields' var]))+abs(T0.(['n_afields' var]))<5,:)=nan;
T1.([crit var '_' crit2])(abs(T1.(['n_fields' var]))+abs(T1.(['n_afields' var]))<5,:)=nan;

end

%%
dat0 = [T0.mRDI_L_UV T0.mRDI_R_UV, T0.mRDI_C_UV];
dat1 = [T1.mRDI_L_UV T1.mRDI_R_UV, T1.mRDI_C_UV];
[m0,t0] = nanmax(abs(dat0),[],2);
[m1,t1] = nanmax(abs(dat1),[],2);

for r=1:size(T0,1)
    T0.mRDI_M_UV(r) = dat0(r,t0(r));
end

for r=1:size(T1,1)
    T1.mRDI_M_UV(r) = dat1(r,t1(r));
end

%% pie_sp and non-sp
crit = 'pBinomDev'; crit2 = 'UV_FR';
R0 = T0(max([T0.nRDI_L_max,T0.nRDI_R_max,T0.nRDI_C_max],[],2)>=5,:);
R1 = T1(max([T1.nRDI_L_max,T1.nRDI_R_max,T1.nRDI_C_max],[],2)>=5,:);

R0 = T0(max([T0.nRDI_L_max,T0.nRDI_R_max,T0.nRDI_C_max],[],2)>=5,:);
% % R0 = R0(R0.rat==232,:);
R1 = T1(max([T1.nRDI_L_max,T1.nRDI_R_max,T1.nRDI_C_max],[],2)>=5,:);

x0 = size(T0,1)-size(R0,1);
x1 = sum(min([R0.([crit '_L_' crit2]),R0.([crit '_R_' crit2]),R0.([crit '_C_' crit2])],[],2)<al & ~(R0.DecodingP_all<al)); 
x2 = sum((min([R0.([crit '_L_' crit2]),R0.([crit '_R_' crit2]),R0.([crit '_C_' crit2])],[],2)<al) & (R0.DecodingP_all<al));
x3 = sum(~(min([R0.([crit '_L_' crit2]),R0.([crit '_R_' crit2]),R0.([crit '_C_' crit2])],[],2)<al) & (R0.DecodingP_all<al)); 
x4 = sum(~(min([R0.([crit '_L_' crit2]),R0.([crit '_R_' crit2]),R0.([crit '_C_' crit2])],[],2)<al) & ~(R0.DecodingP_all<al));
dat0 = [x1 x2 x3 x4];

y0 = size(T1,1)-size(R1,1);
y1 = sum(min([R1.([crit '_L_' crit2]),R1.([crit '_R_' crit2]),R1.([crit '_C_' crit2])],[],2)<al & ~(R1.DecodingP_all<al)); 
y2 = sum((min([R1.([crit '_L_' crit2]),R1.([crit '_R_' crit2]),R1.([crit '_C_' crit2])],[],2)<al) & (R1.DecodingP_all<al));
y3 = sum(~(min([R1.([crit '_L_' crit2]),R1.([crit '_R_' crit2]),R1.([crit '_C_' crit2])],[],2)<al) & (R1.DecodingP_all<al)); 
y4 = sum(~(min([R1.([crit '_L_' crit2]),R1.([crit '_R_' crit2]),R1.([crit '_C_' crit2])],[],2)<al) & ~(R1.DecodingP_all<al)); 
dat1 = [y1 y2 y3 y4];

figure
subplot(1,2,1)
pie(dat0); title(['SUB - ' num2str(x4) ',' num2str(x3) ',' num2str(x2) ',' num2str(x1) ',' num2str(x0)])
subplot(1,2,2)
pie(dat1); title(['CA1 - ' num2str(y4) ',' num2str(y3) ',' num2str(y2) ',' num2str(y1) ',' num2str(y0)])


%% pie_scene and Choice
crit = 'pBinomDev'; crit2 = 'UV_FR';
R0 = T0(max([T0.nRDI_L_max,T0.nRDI_R_max,T0.nRDI_C_max],[],2)>=5,:);
% % R0 = R0(R0.rat==232,:);
R1 = T1(max([T1.nRDI_L_max,T1.nRDI_R_max,T1.nRDI_C_max],[],2)>=5,:);
% R1 = R1(R1.session~=1,:);

% R0 = T0(T0.nFields>=5,:);
% R1 = T1(T1.nFields>=5,:);

x1 = sum((min([R0.([crit '_L_' crit2]),R0.([crit '_R_' crit2])],[],2)<al) & ~(min([R0.([crit '_C_' crit2])],[],2)<al)); 
x2 = sum((min([R0.([crit '_L_' crit2]),R0.([crit '_R_' crit2])],[],2)<al) & (min([R0.([crit '_C_' crit2])],[],2)<al)); 
x3 = sum(~(min([R0.([crit '_L_' crit2]),R0.([crit '_R_' crit2])],[],2)<al) & (min([R0.([crit '_C_' crit2])],[],2)<al)); 
x4 = sum(~(min([R0.([crit '_L_' crit2]),R0.([crit '_R_' crit2])],[],2)<al) & ~(min([R0.([crit '_C_' crit2])],[],2)<al)); 

x5 = sum(~(R0.([crit '_L_' crit2])<al) & (R0.([crit '_R_' crit2])<al) & (R0.([crit '_C_' crit2])<al)); 
x6 = sum((R0.([crit '_L_' crit2])<al) & ~(R0.([crit '_R_' crit2])<al) & (R0.([crit '_C_' crit2])<al)); 
x7 = sum((R0.([crit '_L_' crit2])<al) & (R0.([crit '_R_' crit2])<al) & ~(R0.([crit '_C_' crit2])<al)); 
x8 = sum((R0.([crit '_L_' crit2])<al) & (R0.([crit '_R_' crit2])<al) & (R0.([crit '_C_' crit2])<al)); 
dat0 = [x1 x2 x3 x4]; dat01 = [x5 x6 x7 x8];
N0 = [x1 x2 x3];

y1 = sum((min([R1.([crit '_L_' crit2]),R1.([crit '_R_' crit2])],[],2)<al) & ~(min([R1.([crit '_C_' crit2])],[],2)<al)); 
y2 = sum((min([R1.([crit '_L_' crit2]),R1.([crit '_R_' crit2])],[],2)<al) & (min([R1.([crit '_C_' crit2])],[],2)<al)); 
y3 = sum(~(min([R1.([crit '_L_' crit2]),R1.([crit '_R_' crit2])],[],2)<al) & (min([R1.([crit '_C_' crit2])],[],2)<al)); 
y4 = sum(~(min([R1.([crit '_L_' crit2]),R1.([crit '_R_' crit2])],[],2)<al) & ~(min([R1.([crit '_C_' crit2])],[],2)<al)); 

y5 = sum(~(R1.([crit '_L_' crit2])<al) & (R1.([crit '_R_' crit2])<al) & (R1.([crit '_C_' crit2])<al)); 
y6 = sum((R1.([crit '_L_' crit2])<al) & ~(R1.([crit '_R_' crit2])<al) & (R1.([crit '_C_' crit2])<al)); 
y7 = sum((R1.([crit '_L_' crit2])<al) & (R1.([crit '_R_' crit2])<al) & ~(R1.([crit '_C_' crit2])<al)); 
y8 = sum((R1.([crit '_L_' crit2])<al) & (R1.([crit '_R_' crit2])<al) & (R1.([crit '_C_' crit2])<al)); 
dat1 = [y1 y2 y3 y4]; dat11 = [y5 y6 y7 y8];
N1 = [y1 y2 y3];

figure
subplot(3,4,1)
pie(dat0); title(['SUB - ' num2str(x4) ',' num2str(x3) ',' num2str(x2) ',' num2str(x1)])
subplot(3,4,2)
pie(N0); 
subplot(3,4,3)
pie(dat1); title(['CA1 - ' num2str(y4) ',' num2str(y3) ',' num2str(y2) ',' num2str(y1)])
subplot(3,4,4)
pie(N1); 

subplot(3,2,3)
bar([[x1 x2 x3]./sum(N0)]); title(['SUB'])

subplot(3,2,4)
bar([[y1 y2 y3]./sum(N1)]); title(['CA1'])

subplot(3,2,5)
bar(dat01 ./sum(N0)); title(['SUB'])
ylim([0 0.25])

subplot(3,2,6)
bar(dat11 ./sum(N1)); title(['CA1'])
ylim([0 0.25])
% ylim([0 .35])

y51 = sum((R1.([crit '_L_' crit2])<al) & (R1.([crit '_C_' crit2])<al) & (R1.(['mRDI_C_' crit2])>0));
y61 = sum((R1.([crit '_R_' crit2])<al) & (R1.([crit '_C_' crit2])<al) & (R1.(['mRDI_C_' crit2])<0));

y62 = sum((R1.([crit '_R_' crit2])<al) & (R1.([crit '_C_' crit2])<al) & (R1.(['mRDI_C_' crit2])>0));
y52 = sum((R1.([crit '_L_' crit2])<al) & (R1.([crit '_C_' crit2])<al) & (R1.(['mRDI_C_' crit2])<0));

x51 = sum((R0.([crit '_L_' crit2])<al) & (R0.([crit '_C_' crit2])<al) & (R0.(['mRDI_C_' crit2])>0));
x61 = sum((R0.([crit '_R_' crit2])<al) & (R0.([crit '_C_' crit2])<al) & (R0.(['mRDI_C_' crit2])<0));

x52 = sum((R0.([crit '_L_' crit2])<al) & (R0.([crit '_C_' crit2])<al) & (R0.(['mRDI_C_' crit2])<0));
x62 = sum((R0.([crit '_R_' crit2])<al) & (R0.([crit '_C_' crit2])<al) & (R0.(['mRDI_C_' crit2])>0));

figure
bar([[x51+x61, x52+x62]/sum(N0); [y51+y61, y52+y62]/sum(N1)])

%% cdf_R^2
R0=T0; R1=T1;
% R0 = T0(T0.DecodingP_all<al,:);
% R1 = T1(T1.DecodingP_all<al,:);
figure; hold on
c0 = cdfplot((R0.Decoding_Rsq));
c0.Color = CList(1,:);
c1 = cdfplot((R1.Decoding_Rsq));
c1.Color = CList(2,:);
[h,p] = kstest2((R0.Decoding_Rsq),(R1.Decoding_Rsq))

%% histogram_R^2
R0=T0; R1=T1;
% R0 = T0(T0.DecodingP_all<al,:);
% R1 = T1(T1.DecodingP_all<al,:);
figure; hold on
h0 = histogram((R0.Decoding_Rsq),'binwidth',0.02,'Normalization','probability');
h0.EdgeColor= CList(1,:); h0.FaceAlpha= 0;
h1 = histogram((R1.Decoding_Rsq),'binwidth',0.02,'Normalization','probability');
h1.EdgeColor= CList(2,:); h1.FaceAlpha= 0;
[h,p] = kstest2((R0.Decoding_Rsq),(R1.Decoding_Rsq))
xlim([0 1])

%% bar_R^2
R0 = T0(T0.DecodingP_all<al,:);
R1 = T1(T1.DecodingP_all<al,:);
figure; hold on
x0 = R0.Decoding_Rsq;
x1 = R1.Decoding_Rsq;
dat = [nanmean(x0) nanmean(x1)];
err = [nanstd(x0)/sqrt(sum(~isnan(x0))) nanstd(x1)/sqrt(sum(~isnan(x1)))];
bar(dat)
errorbar(dat,err)
% yticks([0:al:0.2])
[p,h,stats] = ranksum(x0,x1)
title(['p=' jjnum2str(p,3) ', t=' jjnum2str(stats.zval,3)])

%% scatter_R^2 and p-value
figure; hold on
x0 = R0.Decoding_Rsq; y0 = R0.DecodingP_all;
x1 = R1.Decoding_Rsq;  y1 = R1.DecodingP_all;
scatter(x0,y0,40,CList(1,:),'filled')
scatter(x1,y1,40,CList(2,:),'filled')
xlabel('R^2 for position decoding')
ylabel('p-value for permutation test')
legend({'SUB','CA1'})
[h,p,~,stats] = ttest2(x0,x1);
% title(['p=' jjnum2str(p,3) ', t=' jjnum2str(stats.tstat,3)])
%%
x1 = sum(T0.mRDI_M_UV==0); x2 = sum(min([T0.pRDI_L_UV,T0.pRDI_R_UV,T0.pRDI_C_UV],[],2)<al & T0.DecodingP_all>=al);
dat0 = [x2, size(T0,1)-x1-x2,x1];
y1 = sum(T1.mRDI_M_UV==0); y2 = sum(min([T1.pRDI_L_UV,T1.pRDI_R_UV,T1.pRDI_C_UV],[],2)<al & T1.DecodingP_all>=al);
dat1 = [ y2, size(T1,1)-y1-y2,y1];

figure
subplot(1,2,1)
pie(dat0); title(['SUB - ' num2str(dat0(3)) ',' num2str(dat0(2)) ',' num2str(dat0(1))])
subplot(1,2,2)
pie(dat1); title(['CA1 - ' num2str(dat1(3)) ',' num2str(dat1(2)) ',' num2str(dat1(1))])

%%
figure;

R0 = T0(min([T0.pRDI_L_UV,T0.pRDI_R_UV],[],2)<al,:);
R1 = T1(min([T1.pRDI_L_UV,T1.pRDI_R_UV],[],2)<al,:);
subplot(2,2,1)
pie([size(T0,1)-size(R0,1),size(R0,1)],'%.3f%%')
title([num2str(size(R0,1)) ', ' num2str(size(T0,1) - size(R0,1))])
subplot(2,2,2)
pie([size(T1,1)-size(R1,1),size(R1,1)],'%.3f%%')
title([num2str(size(R1,1)) ', ' num2str(size(T1,1) - size(R1,1))])

R0 = T0(min([T0.pRDI_C_UV],[],2)<al,:);
R1 = T1(min([T1.pRDI_C_UV],[],2)<al,:);
subplot(2,2,3)
pie([size(T0,1)-size(R0,1),size(R0,1)],'%.3f%%')
title([num2str(size(R0,1)) ', ' num2str(size(T0,1) - size(R0,1))])
subplot(2,2,4)
pie([size(T1,1)-size(R1,1),size(R1,1)],'%.3f%%')
title([num2str(size(R1,1)) ', ' num2str(size(T1,1) - size(R1,1))])

%% 
crit = 'pBinomDev'; crit2 = 'UV';
varList  = {'_L_' ,'_R_' , '_C_' };
figure;
for v=1:3
    %% pie_L/R/C SI
var = varList{v};
% R0 = T0(T0.(['p' crit var])<al,:);
% R1 = T1(T1.(['p' crit var])<al,:);

% R0 = T0(abs(T0.(['n_fields' var]))+abs(T0.(['n_afields' var]))>=5,:);
% R1 = T1(abs(T1.(['n_fields' var]))+abs(T1.(['n_afields' var]))>=5,:);

% R0 = T0; R1=T1;


R0 = T0(T0.(['nRDI' var 'max'])>=5,:);
% R0 = R0(R0.rat==232,:);
R1 = T1(T1.(['nRDI' var 'max'])>=5,:);

S0 = R0(R0.([crit var crit2])<al,:);
S1 = R1(R1.([crit var crit2])<al,:);


subplot(4,6,v*2-1)
pie([size(R0,1)-size(S0,1),size(S0,1)],'%.3f%%')
title([num2str(size(S0,1)) ', ' num2str(size(R0,1)-size(S0,1))])
subplot(4,6,v*2)
pie([size(R1,1)-size(S1,1),size(S1,1)],'%.3f%%')
title([num2str(size(S1,1)) ', '  num2str(size(R1,1)-size(S1,1))])

%% cdf_L/R/C SI
subplot(4,3,v+3); hold on
c0 = cdfplot(abs(R0.(['mRDI' var crit2])));
c0.Color = CList(1,:);
c1 = cdfplot(abs(R1.(['mRDI' var crit2])));
c1.Color = CList(2,:);
[h,p,stats] = kstest2((R0.(['mRDI' var crit2])),(R1.(['mRDI' var crit2])))
title(['p=' jjnum2str(p,4) ',k=' jjnum2str(stats,4)])
xlim([0 0.8])

%% bar_L/R/C SI
subplot(4,3,v+6); hold on
x0 = abs(S0.(['mRDI' var crit2]));
x1 = abs(S1.(['mRDI' var crit2]));
dat = [nanmean(x0) nanmean(x1)];
err = [nanstd(x0)/sqrt(sum(~isnan(x0))) nanstd(x1)/sqrt(sum(~isnan(x1)))];
bar(dat)
errorbar(dat,err)
% yticks([0:al:0.2])
[p,h,stats] = ranksum(x0,x1)
title(['p=' jjnum2str(p,3) ', Z=' jjnum2str(stats.zval,3)])

%% histogram_L/R/C SI

subplot(4,3,v+9); hold on



h0 = histogram(abs(S0.(['mRDI' var  crit2])),'binwidth',0.02,'Normalization','probability');
h0.EdgeColor= CList(1,:); h0.FaceAlpha= 0;
h1 = histogram(abs(S1.(['mRDI' var  crit2])),'binwidth',0.02,'Normalization','probability');
h1.EdgeColor= CList(2,:); h1.FaceAlpha= 0;


% [h,p] = kstest2((R0.(['m' var])),(R1.(['m' var])))
xlim([0 0.8])


% Fit to a normal (or a different distribution if you choose)
pd = fitdist(abs(S0.(['mRDI' var  crit2])),'Kernel');
% Find the pdf that spans the disribution
x_pdf = linspace(0,.8);
y_pdf = pdf(pd,x_pdf);
% Plot
line(x_pdf,y_pdf/(40*0.8),'LineWidth',2,'color',CList(1,:));

% Fit to a normal (or a different distribution if you choose)
pd = fitdist(abs(S1.(['mRDI' var  crit2])),'Kernel');
% Find the pdf that spans the disribution
x_pdf = linspace(0,.8);
y_pdf = pdf(pd,x_pdf);
% Plot
line(x_pdf,y_pdf/sum(y_pdf),'LineWidth',2,'color',CList(2,:));

end
%% cdf_max SI
% R0 = T0(min([T0.pRDI_L_UV,T0.pRDI_R_UV,T0.pRDI_C_UV],[],2)<al,:);
% R1 = T1(min([T1.pRDI_L_UV,T1.pRDI_R_UV,T1.pRDI_C_UV],[],2)<al,:);
figure; hold on
c0 = cdfplot(abs(R0.mRDI_M_UV));
c0.Color = CList(1,:);
c1 = cdfplot(abs(R1.mRDI_M_UV));
c1.Color = CList(2,:);
[h,p] = kstest2((R0.mRDI_M_UV),(R1.mRDI_M_UV))
%% boxplot_max SI
figure
R0 = T0(min([T0.pRDI_L_UV,T0.pRDI_R_UV,T0.pRDI_C_UV],[],2)<al,:);
R1 = T1(min([T1.pRDI_L_UV,T1.pRDI_R_UV,T1.pRDI_C_UV],[],2)<al,:);
    x1 = abs([R0.mRDI_M_UV]);
    x2 = abs([R1.mRDI_M_UV]);
x = [x1;x2];
g = [repmat({['SUB (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['CA1 (n=' num2str(size(x2,1)) ')']},size(x2,1),1)];
boxplot(x,g)
ylabel('ensemble mean of selectivity index')
xlabel('max')
% ylim([0 1])
[~,p] = ttest2(x1,x2)
%% bar_max SI
R0 = T0(min([T0.pRDI_L_UV,T0.pRDI_R_UV,T0.pRDI_C_UV],[],2)<al,:);
R1 = T1(min([T1.pRDI_L_UV,T1.pRDI_R_UV,T1.pRDI_C_UV],[],2)<al,:);
figure; hold on
x0 = abs(R0.mRDI_M_UV);
x1 = abs(R1.mRDI_M_UV);
dat = [nanmean(x0) nanmean(x1)];
err = [nanstd(x0)/sqrt(sum(~isnan(x0))) nanstd(x1)/sqrt(sum(~isnan(x1)))];
bar(dat)
errorbar(dat,err)
% yticks([0:al:0.2])
[p,h,stats] = ranksum(x0,x1)
title(['p=' jjnum2str(p,3) ', t=' jjnum2str(stats.zval,3)])

%% bar_Ratio
crit = 'pRatio'; crit2 = 'UV';
R0 = T0(max([T0.nRDI_L_max,T0.nRDI_R_max,T0.nRDI_C_max],[],2)>=5,:);
R1 = T1(max([T1.nRDI_L_max,T1.nRDI_R_max,T1.nRDI_C_max],[],2)>=5,:);


R0.([crit '_M_' crit2]) = max([R0.([crit '_L_' crit2]),R0.([crit '_R_' crit2]),R0.([crit '_C_' crit2])],[],2);
R1.([crit '_M_' crit2]) = max([R1.([crit '_L_' crit2]),R1.([crit '_R_' crit2]),R1.([crit '_C_' crit2])],[],2);

figure; hold on
x0 = abs(R0.([crit '_M_' crit2]));
x1 = abs(R1.([crit '_M_' crit2]));
dat = [nanmean(x0) nanmean(x1)];
err = [nanstd(x0)/sqrt(sum(~isnan(x0))) nanstd(x1)/sqrt(sum(~isnan(x1)))];
bar(dat)
errorbar(dat,err)
% yticks([0:al:0.2])
[p,h,stats] = ranksum(x0,x1)
title(['p=' jjnum2str(p,3) ', t=' jjnum2str(stats.zval,3)])

%% bar_Ratio_nsp non_nsp
crit = 'pRatio'; crit2 = 'UV';
R0 = T0(max([T0.nRDI_L_max,T0.nRDI_R_max,T0.nRDI_C_max],[],2)>=5,:);
R1 = T1(max([T1.nRDI_L_max,T1.nRDI_R_max,T1.nRDI_C_max],[],2)>=5,:);

R0.([crit '_M_' crit2]) = max([R0.([crit '_L_' crit2]),R0.([crit '_R_' crit2]),R0.([crit '_C_' crit2])],[],2);
R1.([crit '_M_' crit2]) = max([R1.([crit '_L_' crit2]),R1.([crit '_R_' crit2]),R1.([crit '_C_' crit2])],[],2);

crit3 = 'pBinomDev'; crit2 = 'UV';
id0 = min([R0.([crit3 '_L_' crit2]),R0.([crit3 '_R_' crit2]),R0.([crit3 '_C_' crit2])],[],2)<al;
id1 = min([R1.([crit3 '_L_' crit2]),R1.([crit3 '_R_' crit2]),R1.([crit3 '_C_' crit2])],[],2)<al;

figure; hold on
x0 = abs(R0.([crit '_M_' crit2])(id0));
x1 = abs(R1.([crit '_M_' crit2])(id1));

y0 = abs(R0.([crit '_M_' crit2])(~id0));
y1 = abs(R1.([crit '_M_' crit2])(~id1));

subplot(1,2,1)
hold on
dat = [nanmean(x0) nanmean(x1)];
err = [nanstd(x0)/sqrt(sum(~isnan(x0))) nanstd(x1)/sqrt(sum(~isnan(x1)))];
bar(dat)
errorbar(dat,err)
[p,h,stats] = ranksum(x0,y0)
title(['p=' jjnum2str(p,3) ', t=' jjnum2str(stats.zval,3)])

subplot(1,2,2)
hold on
dat = [nanmean(y0) nanmean(y1)];
err = [nanstd(y0)/sqrt(sum(~isnan(y0))) nanstd(y1)/sqrt(sum(~isnan(y1)))];
bar(dat)
errorbar(dat,err)
ylim([0 1])
% yticks([0:al:0.2])
[p,h,stats] = ranksum(y1,x1)
title(['p=' jjnum2str(p,3) ', t=' jjnum2str(stats.zval,3)])
%% scatter_MSI and p-value
figure; hold on
R0 = T0(min([T0.pRDI_L_UV,T0.pRDI_R_UV,T0.pRDI_C_UV],[],2)<al,:);
R1 = T1(min([T1.pRDI_L_UV,T1.pRDI_R_UV,T1.pRDI_C_UV],[],2)<al,:);
x0 = abs(R0.mRDI_M_UV); y0 = R0.DecodingP_all;
x1 = abs(R1.mRDI_M_UV);  y1 = R1.DecodingP_all;
scatter(x0,y0,40,CList(1,:),'filled')
scatter(x1,y1,40,CList(2,:),'filled')
xlabel('Abs. ensemble mean of selectivity index (max)')
ylabel('p-value for permutation test')
legend({'SUB','CA1'})
%%
figure; hold on
x = [1 2 3 4];
m0 = nanmax(abs([T0.mRDI_L_UV T0.mRDI_R_UV, T0.mRDI_C_UV]),[],2);
m1 = nanmax(abs([T1.mRDI_L_UV T1.mRDI_R_UV, T1.mRDI_C_UV]),[],2);
data = [nanmean(abs(T0.mRDI_L_UV)) nanmean(abs(T1.mRDI_L_UV));...
    nanmean(abs(T0.mRDI_R_UV)) nanmean(abs(T1.mRDI_R_UV));...
    nanmean(abs(T0.mRDI_C_UV)) nanmean(abs(T1.mRDI_C_UV));...
    nanmean(abs(m0)) nanmean(abs(m1))];
err = [nanstd(abs(T0.mRDI_L_UV))/sqrt(sum(~isnan(T0.mRDI_L_UV))) nanstd(abs(T1.mRDI_L_UV))/sqrt(sum(~isnan(T1.mRDI_L_UV)));...
    nanstd(abs(T0.mRDI_R_UV))/sqrt(sum(~isnan(T0.mRDI_R_UV))) nanstd(abs(T1.mRDI_R_UV))/sqrt(sum(~isnan(T1.mRDI_R_UV)));...
    nanstd(abs(T0.mRDI_C_UV))/sqrt(sum(~isnan(T0.mRDI_C_UV))) nanstd(abs(T1.mRDI_C_UV))/sqrt(sum(~isnan(T1.mRDI_C_UV)));...
    nanstd(abs(m0))/sqrt(sum(~isnan(m0))) nanstd(abs(m1))/sqrt(sum(~isnan(m1)))];

b = bar(x,data);
b(1).FaceColor = CList(1,:);
b(2).FaceColor = CList(2,:);

for k = 1:size(data,2)
    % get x positions per group
    xpos = b(k).XData + b(k).XOffset;
    % draw errorbar
    errorbar(xpos, data(:,k), err(:,k), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
end
legend({'SUB','CA1'},'location','eastoutside')
set(gca,'fontsize',12,'fontweight','b')
ylabel('Ensemble mean of selectivity index')
xticks([1 2 3 4])
xticklabels({'L Scene', 'R Scene','Choice', 'max.'})
xlim([3.5 4.5])

x=abs(abs(m0)); y=abs(abs(m1));
x=abs(T0.mRDI_C_UV); y = abs(T1.mRDI_C_UV);
[p,h,~,stats] = ttest2(x,y)

%% boxplot
figure
%     x1 = nanmax(abs([T0.mRDI_L_UV T0.mRDI_R_UV, T0.mRDI_C_UV]),[],2);
%     x2 = nanmax(abs([T1.mRDI_L_UV T1.mRDI_R_UV, T1.mRDI_C_UV]),[],2);

    x1 = nanmax(abs([T0.nRDI_L_UV]),[],2);
    x2 = nanmax(abs([T1.nRDI_L_UV]),[],2);
x = [x1;x2];
g = [repmat({['SUB (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['CA1 (n=' num2str(size(x2,1)) ')']},size(x2,1),1)];
boxplot(x,g)
ylabel('ensemble mean of selectivity index')
xlabel('left')
% ylim([0 1])
[~,p] = ttest2(x1,x2)

%%
figure; hold on
%     x1 = nanmax(abs([T0.mRDI_L_UV T0.mRDI_R_UV, T0.mRDI_C_UV]),[],2);
%     x2 = nanmax(abs([T1.mRDI_L_UV T1.mRDI_R_UV, T1.mRDI_C_UV]),[],2);

    x1 = nanmax(abs([T0.nRDI_C_UV]),[],2);
    x2 = nanmax(abs([T1.nRDI_C_UV]),[],2);
histogram(x1,'facecolor',CList(1,:))
histogram(x2,'facecolor',CList(2,:))
xlabel('Rate remapping fields (L Scene)')
ylabel('# of ripples')
% ylim([0 1])
[~,p] = ttest2(x1,x2)
legend({'SUB (n=1746)','CA1 (n=1534)'})
%%
m0 = nanmin(abs([T0.nRDI_L_UV T0.nRDI_R_UV, T0.nRDI_C_UV]),[],2);
m1 = nanmin(abs([T1.nRDI_L_UV T1.nRDI_R_UV, T1.nRDI_C_UV]),[],2);
sum(m1<5)
%%
T0.pRDI_L_UV(isnan(T0.pRDI_L_UV))=01;
T1.pRDI_L_UV(isnan(T1.pRDI_L_UV))=01;
T0.pRDI_R_UV(isnan(T0.pRDI_R_UV))=01;
T1.pRDI_R_UV(isnan(T1.pRDI_R_UV))=01;
T0.pRDI_C_UV(isnan(T0.pRDI_C_UV))=1;
T1.pRDI_C_UV(isnan(T1.pRDI_C_UV))=1;
%%
figure; hold on
x = [1 2 3 4];
m0 = nanmin([T0.pRDI_L_UV T0.pRDI_R_UV, T0.pRDI_C_UV],[],2);
m1 = nanmin([T1.pRDI_L_UV T1.pRDI_R_UV, T1.pRDI_C_UV],[],2);
data = [nanmean(abs(T0.pRDI_L_UV)) nanmean(abs(T1.pRDI_L_UV));...
    nanmean(abs(T0.pRDI_R_UV)) nanmean(abs(T1.pRDI_R_UV));...
    nanmean(abs(T0.pRDI_C_UV)) nanmean(abs(T1.pRDI_C_UV));
    nanmean(m0) nanmean(m1)];
err = [nanstd(abs(T0.pRDI_L_UV))/sqrt(sum(~isnan(T0.pRDI_L_UV))) nanstd(abs(T1.pRDI_L_UV))/sqrt(sum(~isnan(T1.pRDI_L_UV)));...
    nanstd(abs(T0.pRDI_R_UV))/sqrt(sum(~isnan(T0.pRDI_R_UV))) nanstd(abs(T1.pRDI_R_UV))/sqrt(sum(~isnan(T1.pRDI_R_UV)));...
    nanstd(abs(T0.pRDI_C_UV))/sqrt(sum(~isnan(T0.pRDI_C_UV))) nanstd(abs(T1.pRDI_C_UV))/sqrt(sum(~isnan(T1.pRDI_C_UV)));
    nanstd(abs(m0))/sqrt(sum(~isnan(m0))) nanstd(abs(m1))/sqrt(sum(~isnan(m1)))];

b = bar(x,data);
b(1).FaceColor = CList(1,:);
b(2).FaceColor = CList(2,:);

for k = 1:size(data,2)
    % get x positions per group
    xpos = b(k).XData + b(k).XOffset;
    % draw errorbar
    errorbar(xpos, data(:,k), err(:,k), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
end
legend({'SUB','CA1'},'location','eastoutside')
set(gca,'fontsize',12,'fontweight','b')
ylabel('p-value')
xticks([1 2 3 4])
xticklabels({'L Scene', 'R Scene','Choice', 'min.'})

x=abs(T0.pRDI_R_UV); y=abs(T1.pRDI_R_UV);
[p,h] = ranksum(x, y)
%%


figure; hold on
x = [1 2];
data = [nanmean(Rip_SUB.NS_p) nanmean(Rip_CA1.NS_p); nanmean(Rip_SUB.S_p) nanmean(Rip_CA1.S_p)];
err = [nanstd(Rip_SUB.NS_p)/sqrt(size(Rip_SUB,1)) nanstd(Rip_CA1.NS_p)/sqrt(size(Rip_CA1,1));...
    nanstd(Rip_SUB.S_p)/sqrt(size(Rip_SUB,1)) nanstd(Rip_CA1.S_p)/sqrt(size(Rip_CA1,1))];

b = bar(x,data);
b(1).FaceColor = CList(1,:);
b(2).FaceColor = CList(2,:);

for k = 1:size(data,1)
    % get x positions per group
    xpos = b(k).XData + b(k).XOffset;
    % draw errorbar
    errorbar(xpos, data(:,k), err(:,k), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
end
legend({'SUB','CA1'})
set(gca,'fontsize',12,'fontweight','b')
ylabel('p-value')
xticks([1 2])
xticklabels({'Non-Spatial','Spatial'})

    %% pie_replay direction
% var = DecodingP_all;
R0 = T0(T0.DecodingP_all<al,:);
R1 = T1(T1.DecodingP_all<al,:);

R0F = R0(R0.DecodingDir>0,:); R0R = R0(R0.DecodingDir<0,:);
R1F = R1(R1.DecodingDir>0,:); R1R = R1(R1.DecodingDir<0,:);

figure;
subplot(2,2,1)
pie([size(T0,1)-size(R0,1),size(R0R,1),size(R0F,1)],'%.3f%%')
title([num2str(size(R0F,1)) ', ' num2str(size(R0R,1)) ', ' num2str(size(T0,1) - size(R0,1))])
subplot(2,2,2)
pie([size(T1,1)-size(R1,1),size(R1R,1),size(R1F,1)],'%.3f%%')
title([num2str(size(R1F,1)) ', ' num2str(size(R1R,1)) ', ' num2str(size(T1,1) - size(R1,1))])

subplot(2,2,3)
bar([[size(R0F,1)]/size(R0,1) [size(R1F,1)]/size(R1,1)]); title(['Forward replay'])
% ylim([0 .2])
subplot(2,2,4)
bar([[size(R0R,1)]/size(R0,1) [size(R1R,1)]/size(R1,1)]); title(['Reverse replay'])
% ylim([0 .1])
%% bar_unittype
fieldName = {'L','R','C'}; varName  = {'Selectivity_LScene','Selectivity_RScene','Selectivity_LR'};

U0 = readtable([ROOT.Units '\UnitsTable_SUB_forAnalysis.xlsx']); U0(U0.rat==232 & U0.session==4,:)=[];
U1 = readtable([ROOT.Units '\UnitsTable_CA1_forAnalysis.xlsx']); U1(U1.rat==232 & U1.session==4,:)=[];

for i=1:3
V0=[]; V1=[]; C0=[]; C1=[];

V0.(fieldName{i})(:,1) = sum((abs(U0.(varName{i}))<1));
V0.(fieldName{i})(:,2) = sum((abs(U0.(varName{i}))==1));
V0.(fieldName{i})(:,3) = sum((abs(U0.(varName{i}))==2));
V0.(fieldName{i})(:,4) = sum((abs(U0.(varName{i}))>2));
% V0.(fieldName{i})(:,4) = sum((abs(U0.(varName{i}))==3));
% V0.(fieldName{i})(:,5) = sum((abs(U0.(varName{i}))==4));


V1.(fieldName{i})(:,1) = sum((abs(U1.(varName{i}))<1));
V1.(fieldName{i})(:,2) = sum((abs(U1.(varName{i}))==1));
V1.(fieldName{i})(:,3) = sum((abs(U1.(varName{i}))==2));
V1.(fieldName{i})(:,4) = sum((abs(U1.(varName{i}))>2));
% V1.(fieldName{i})(:,4) = sum((abs(U1.(varName{i}))==3));
% V1.(fieldName{i})(:,5) = sum((abs(U1.(varName{i}))==4));

C0.(fieldName{i})(:,1) = (T0.nPCs-sum([T0.(['n_sf_' fieldName{i}]) ,T0.(['n_mfs_' fieldName{i}]) ,T0.(['n_hom_' fieldName{i}]) ,T0.(['n_het_' fieldName{i}]) ],2)) ./ T0.nPCs;
C0.(fieldName{i})(:,2) = T0.(['n_sf_' fieldName{i}]) ./ T0.nPCs;
C0.(fieldName{i})(:,3) = T0.(['n_mfs_' fieldName{i}])  ./ T0.nPCs;
C0.(fieldName{i})(:,4) = (T0.(['n_hom_' fieldName{i}]) +T0.(['n_het_' fieldName{i}]) ) ./ T0.nPCs;
% C0.(fieldName{i})(:,4) = T0.(['n_hom_' fieldName{i}]) ./ T0.nPCs;
% C0.(fieldName{i})(:,5) = T0.(['n_het_' fieldName{i}]) ./ T0.nPCs;

C1.(fieldName{i})(:,1) = (T1.nPCs-sum([T1.(['n_sf_' fieldName{i}]) ,T1.(['n_mfs_' fieldName{i}]) ,T1.(['n_hom_' fieldName{i}]) ,T1.(['n_het_' fieldName{i}]) ],2)) ./ T1.nPCs;
C1.(fieldName{i})(:,2) = T1.(['n_sf_' fieldName{i}]) ./ T1.nPCs;
C1.(fieldName{i})(:,3) = T1.(['n_mfs_' fieldName{i}])  ./ T1.nPCs;
C1.(fieldName{i})(:,4) = (T1.(['n_hom_' fieldName{i}]) +T1.(['n_het_' fieldName{i}]) ) ./ T1.nPCs;
% C1.(fieldName{i})(:,4) = T1.(['n_hom_' fieldName{i}]) ./ T1.nPCs;
% C1.(fieldName{i})(:,5) = T1.(['n_het_' fieldName{i}]) ./ T1.nPCs;

figure('position',[111,246,1161,642]);
sgtitle(varName{i})
subplot(2,3,1)
boxplot(C0.(fieldName{i})); ylim([0 1])

subplot(2,3,2); hold on
x0 = abs(C0.(fieldName{i}));
dat = [nanmean(x0)];
err = [nanstd(x0)./sqrt(size(x0,1))];
b=bar(dat);
errorbar(dat,err)
b.FaceColor = CList(1,:);
xticks([1:4]); xticklabels({'N','SF','MF,s','MF,m'})
ylim([0 1])
for d=1:4, [h,p,stats] = signrank(x0(:,d)-V0.(fieldName{i})(d)/size(U0,1)); z(d) = stats.zval; hval(i*2-1,d)=h; pval{d} = jjnum2str(h,4); end


subplot(2,3,3)
b=bar(V0.(fieldName{i})/size(U0,1),'stacked'); b.FaceColor = CList(1,:);
ylim([0 1])
title(['SUB, unit' '-' num2str(hval(i*2-1,:))])

subplot(2,3,1)
bar([V0.(fieldName{i})/size(U0,1);dat],'stacked')
title(['SUB, ripple - ' num2str(z) ])
ylim([0 1])

subplot(2,3,4)
boxplot(C1.(fieldName{i})); ylim([0 1])
xticks([1:4]); xticklabels({'N','SF','MF,s','MF,m'})

subplot(2,3,5); hold on
x0 = abs(C1.(fieldName{i}));
dat = [nanmean(x0)];
err = [nanstd(x0)./sqrt(size(x0,1))];
b=bar(dat);
errorbar(dat,err)
b.FaceColor = CList(2,:);

xticks([1:4]); xticklabels({'N','SF','MF,s','MF,m'})
ylim([0 1])
for d=1:4, [h,p,stats] = signrank(x0(:,d)-V1.(fieldName{i})(d)/size(U1,1)); z(d) = stats.zval;  hval(i*2,d)=h; pval{d} = jjnum2str(h,4); end


subplot(2,3,6)
b=bar(V1.(fieldName{i})/size(U1,1)); b.FaceColor = CList(2,:);
xticks([1:4]); xticklabels({'N','SF','MF,s','MF,m'})
ylim([0 1])
title(['CA1, unit' '-' num2str(hval(i*2,:))])


subplot(2,3,4)
bar([V1.(fieldName{i})/size(U1,1);dat],'stacked')
title(['CA1, ripple - ' num2str(z)])
ylim([0 1])

saveas(gca, [ROOT.Processed '/Manuscript/' varName{i} '_celltype.png'])
saveas(gca, [ROOT.Processed '/Manuscript/' varName{i} '_celltype.svg'])
end
%% SWR purity
fieldList  = {'L','R','C'};
varList  = {'RDI_L_' crit2,'RDI_R_' crit2, 'RDI_C_' crit2};
figure;
for v=1:3
    %% pie_L/R/C SI
var = varList{v};
R0 = T0(T0.(['p' var])<al,:);
R1 = T1(T1.(['p' var])<al,:);

thisF = fieldList{v};

subplot(3,6,v*2-1)
pie([size(T0,1)-size(R0,1),size(R0,1)],'%.3f%%')
title([num2str(size(R0,1)) ', ' num2str(size(T0,1) - size(R0,1))])
subplot(3,6,v*2)
pie([size(T1,1)-size(R1,1),size(R1,1)],'%.3f%%')
title([num2str(size(R1,1)) ', ' num2str(size(T1,1) - size(R1,1))])

%% cdf_L/R/C SI
subplot(3,3,v+3); hold on
c0 = cdfplot(abs(R0.(['n_fields_' thisF]))./R0.nFields);
c0.Color = CList(1,:);
c1 = cdfplot(abs(R1.(['n_fields_' thisF]))./R1.nFields);
c1.Color = CList(2,:);
% [h,p] = kstest2((R0.(['m' var])),(R1.(['m' var])))
xlim([0 1])

%% bar_L/R/C SI
subplot(3,9,3*v-2+18); hold on
x0 = abs(R0.(['n_fields_' thisF]))./R0.nFields;
x1 = abs(R1.(['n_fields_' thisF]))./R1.nFields;
dat = [nanmean(x0) nanmean(x1)];
err = [nanstd(x0)/sqrt(sum(~isnan(x0))) nanstd(x1)/sqrt(sum(~isnan(x1)))];
bar(dat)
errorbar(dat,err)
% yticks([0:al:0.2])
[p,h,stats] = ranksum(x0,x1)
title(['p=' jjnum2str(p,3) ', Z=' jjnum2str(stats.zval,3)])
ylim([0 1])

subplot(3,9,3*v-1+18); hold on
x0 = abs(R0.(['n_afields_' thisF]))./R0.nFields;
x1 = abs(R1.(['n_afields_' thisF]))./R1.nFields;
dat = [nanmean(x0) nanmean(x1)];
err = [nanstd(x0)/sqrt(sum(~isnan(x0))) nanstd(x1)/sqrt(sum(~isnan(x1)))];
bar(dat)
errorbar(dat,err)
% yticks([0:al:0.2])
[p,h,stats] = ranksum(x0,x1)
title(['p=' jjnum2str(p,3) ', Z=' jjnum2str(stats.zval,3)])
ylim([0 1])

subplot(3,9,3*v+18); hold on
x0 = 1-(abs(R0.(['n_fields_' thisF]))+abs(R0.(['n_afields_' thisF])))./R0.nFields;
x1 = 1-(abs(R1.(['n_fields_' thisF]))+abs(R1.(['n_afields_' thisF])))./R1.nFields;
dat = [nanmean(x0) nanmean(x1)];
err = [nanstd(x0)/sqrt(sum(~isnan(x0))) nanstd(x1)/sqrt(sum(~isnan(x1)))];
bar(dat)
errorbar(dat,err)
% yticks([0:al:0.2])
[p,h,stats] = ranksum(x0,x1)
title(['p=' jjnum2str(p,3) ', Z=' jjnum2str(stats.zval,3)])
ylim([0 1])


end

%% hist_fields
crit='Binom';
varList  = {'_L','_R', '_C'};
figure;
for v=1:3
    %% pie_L/R/C SI
var = varList{v};
R0 = T0(abs(T0.(['n_fields' var]))+abs(T0.(['n_afields' var]))>=5,:);
R1 = T1(abs(T1.(['n_fields' var]))+abs(T1.(['n_afields' var]))>=5,:);
R0 = R0(R0.(['p' crit var '_' crit2])<al,:);
R1 = R1(R1.(['p' crit var '_' crit2])<al,:);

subplot(3,2,v*2-2+1)
histogram(abs(R0.(['n_fields' var]))+abs(R0.(['n_afields' var])),'binwidth',1,'facecolor',CList(1,:),'Normalization','probability')
xlim([5 17]); ylim([0 .8]); xlabel('number of selective fields'); ylabel('proportion')
title(['n= ' num2str(size(R0,1))])

subplot(3,2,v*2-2+2)
histogram(abs(R1.(['n_fields' var]))+abs(R1.(['n_afields' var])),'binwidth',1,'facecolor',CList(2,:),'Normalization','probability')
xlim([5 17]);  ylim([0 .8]); xlabel('number of selective fields'); ylabel('proportion')
title(['n= ' num2str(size(R1,1))])
end

%% histogram_pBinom
crit = 'pBinomDev'; crit2 = 'UV';
R0 = T0(max([T0.nRDI_L_max,T0.nRDI_R_max,T0.nRDI_C_max],[],2)>=5,:);
% % R0 = R0(R0.rat==232,:);
R1 = T1(max([T1.nRDI_L_max,T1.nRDI_R_max,T1.nRDI_C_max],[],2)>=5,:);

figure
RList={'L','R','C'};
for p=1:3
    rs = RList{p};
    switch p, case 1, xlab='SSI_L'; case 2, xlab='SSI_R'; case 3, xlab='CSI'; end
subplot(2,3,p)
hold on
% c0 = cdfplot(abs(max([abs(R0.mRDI_L_max),abs(R0.mRDI_R_max)],[],2)-abs(R0.mRDI_C_max)));
% c0.Color = CList(1,:);
% c1 = cdfplot(abs(max([abs(R1.mRDI_L_max),abs(R1.mRDI_R_max)],[],2)-abs(R1.mRDI_C_max)));
% c1.Color = CList(2,:);
x1 = R0.([crit '_' rs  '_' crit2]); 
x2 = R1.([crit '_' rs  '_' crit2]); 
% x1 = R0.([crit '_' rs  '_' crit2])(R0.([crit '_' rs  '_' crit2])<0.05); 
% x2 = R0.([crit '_' rs  '_' crit2])(R0.([crit '_' rs  '_' crit2])>=0.05);  
% y1 = R1.([crit '_' rs  '_' crit2])(R1.([crit '_' rs  '_' crit2])<0.05);  
% y2 = R1.([crit '_' rs  '_' crit2])(R1.([crit '_' rs  '_' crit2])>=0.05);  
% [PDF, xValues] = ksdensity(x1); plot(xValues, PDF, 'color',CList(1,:), 'LineWidth', 2);
% [PDF, xValues] = ksdensity(x2); plot(xValues, PDF,'color', CList(1,:), 'linestyle',':','LineWidth', 2);
% [PDF, xValues] = ksdensity(y1); plot(xValues, PDF, 'color',CList(2,:), 'LineWidth', 2);
%  [PDF, xValues] = ksdensity(y2); plot(xValues, PDF,'color', CList(2,:), 'linestyle',':','LineWidth', 2);

histogram(x1,'Binedges',[0:0.02:1],'Facecolor',CList(1,:),'normalization','probability')
histogram(x2,'Binedges',[0:0.02:1],'Facecolor',CList(2,:),'normalization','probability')
xlabel(['p(' xlab ')']); ylim([0 .25]); ylabel('Proportion'); 


subplot(2,3,p+3)
hold on
c0=cdfplot(x1); c0.Color = CList(1,:);
c1=cdfplot(x2); c1.Color = CList(2,:);
xlabel(['p(' xlab ')']); xlim([0 1]); ylabel('Cumulative proportion'); 
end

%% Generate some normally distributed random data
data = randn(1000,1);

% Create a histogram with PDF normalization
figure;
histogram(data, 'Normalization', 'pdf');

% Add labels
xlabel('Data');
ylabel('Probability Density');

% Add a title
title('Normalized PDF Plot of Data');
[PDF, xValues] = ksdensity(data);
hold on;
plot(xValues, PDF, 'r', 'LineWidth', 2);
hold off;

%% scatter_p_Binom
crit = 'pBinomDev'; crit2 = 'UV';
R0 = T0(max([T0.nRDI_L_max,T0.nRDI_R_max,T0.nRDI_C_max],[],2)>=5,:);
% % R0 = R0(R0.rat==232,:);
R1 = T1(max([T1.nRDI_L_max,T1.nRDI_R_max,T1.nRDI_C_max],[],2)>=5,:);

figure
hold on
% c0 = cdfplot(abs(max([abs(R0.mRDI_L_max),abs(R0.mRDI_R_max)],[],2)-abs(R0.mRDI_C_max)));
% c0.Color = CList(1,:);
% c1 = cdfplot(abs(max([abs(R1.mRDI_L_max),abs(R1.mRDI_R_max)],[],2)-abs(R1.mRDI_C_max)));
% c1.Color = CList(2,:);
x1 = abs(min([abs(R0.([crit '_L_' crit2])),abs(R0.([crit '_R_' crit2]))],[],2)); y1=abs(R0.([crit '_C_' crit2]));
x2 = abs(min([abs(R1.([crit '_L_' crit2])),abs(R1.([crit '_R_' crit2]))],[],2)); y2 = abs(R1.([crit '_C_' crit2]));
x1=1-x1;x2=1-x2;y1=1-y1;y2=1-y2;
c0 = cdfplot(abs(x1-y1));
c0.Color = CList(1,:);
c1 = cdfplot(abs(x2-y2));
c1.Color = CList(2,:);
[h,p,stats] = kstest2(x1,x2)
[h,p,stats] = ranksum(x1,x2)
xlim([0 0.8])

%%
figure; hold on
x0 = abs(x1-y1);
x1 = abs(x2-y2);
dat = [nanmean(x0) nanmean(x1)];
err = [nanstd(x0)/sqrt(sum(~isnan(x0))) nanstd(x1)/sqrt(sum(~isnan(x1)))];
bar(dat)
errorbar(dat,err)
% yticks([0:al:0.2])
[p,h,stats] = ranksum(x0,x1)
title(['p=' jjnum2str(p,3) ', t=' jjnum2str(stats.zval,3)])
%%
% x=[x1;x2]; y=[y1;y2];
% g = [zeros(size(x1,1),1);ones(size(x2,1),1)];
% 
% ESD=table;
% ESD.SSI=x; ESD.CSI=y; ESD.region=g;
% [h,a,c,s] = aoctool(x,y,g,0.05,'p(SSI)','p(CSI)',{'SUB','CA1'},'off','separate lines')
% 
figure
scatter(x1,y1,40,CList(1,:),'filled');
x1(isnan(y1))=[]; y1(isnan(y1))=[]; y1(isnan(x1))=[]; x1(isnan(x1))=[];
corrcoef(x1,y1)
hold on
scatter(x2,y2,40,CList(2,:),'filled');
x2(isnan(y2))=[]; y2(isnan(y2))=[]; y2(isnan(x2))=[]; x2(isnan(x2))=[];
corrcoef(x2,y2)
xlim([0 1])

xlabel('p(SSI)'); ylabel('p(CSI)')
set(gca,'fontsize',12','fontweight','b')

figure
hold on
histogram(x1)
histogram(y1)


%% scatter_p_SP and p_NonSP

crit = 'pBinomDev'; crit2 = 'UV';
R0 = T0(max([T0.nRDI_L_max,T0.nRDI_R_max,T0.nRDI_C_max],[],2)>=5,:);
% % R0 = R0(R0.rat==232,:);
R1 = T1(max([T1.nRDI_L_max,T1.nRDI_R_max,T1.nRDI_C_max],[],2)>=5,:);

% c0 = cdfplot(abs(max([abs(R0.mRDI_L_max),abs(R0.mRDI_R_max)],[],2)-abs(R0.mRDI_C_max)));
% c0.Color = CList(1,:);
% c1 = cdfplot(abs(max([abs(R1.mRDI_L_max),abs(R1.mRDI_R_max)],[],2)-abs(R1.mRDI_C_max)));
% c1.Color = CList(2,:);
x1 = abs(min([abs(R0.([crit '_L_' crit2])),abs(R0.([crit '_R_' crit2])) abs(R0.([crit '_C_' crit2]))],[],2));
y1 = R0.DecodingP_all;
x2 = abs(min([abs(R1.([crit '_L_' crit2])),abs(R1.([crit '_R_' crit2])) abs(R1.([crit '_C_' crit2]))],[],2)); 
y2 = R1.DecodingP_all;

figure;
scatter(x1,y1,40,CList(1,:),'filled');
x1(isnan(y1))=[]; y1(isnan(y1))=[]; y1(isnan(x1))=[]; x1(isnan(x1))=[];
corrcoef(x1,y1)
hold on
scatter(x2,y2,40,CList(2,:),'filled');
x2(isnan(y2))=[]; y2(isnan(y2))=[]; y2(isnan(x2))=[]; x2(isnan(x2))=[];
corrcoef(x2,y2)
xlim([0 1]); ylim([0 1])

xlabel('p(SI)'); ylabel('p(Position)')
set(gca,'fontsize',12','fontweight','b')


%%
% Generate bivariate normal distributed data
mu = [0 0]; % mean
sigma = [1 0.5; 0.5 2]; % covariance matrix
data = mvnrnd(mu, sigma, 1000); % generate 1000 data points

% Create scatter plot
figure;
scatter(data(:,1), data(:,2), 'b.');

% Compute the eigenvalues and eigenvectors of the covariance matrix
[V,D] = eig(sigma);

% Compute the angle of tilt of the ellipse
theta = atan2(V(2,1), V(1,1));

% Compute the semi-axes of the ellipse
a = sqrt(5.991 * D(1,1)); % 5.991 corresponds to 95% coverage for 2 variables
b = sqrt(5.991 * D(2,2)); % 5.991 corresponds to 95% coverage for 2 variables

% Draw the 95% coverage region (ellipse)
hold on;
scircle1(a, b, theta, mu(1), mu(2), 'r');
hold off;