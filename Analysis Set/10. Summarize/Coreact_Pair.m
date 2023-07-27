Initial_SWRFilter_common;
warning off

ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip5 = [ROOT.Mother '\Processed Data\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R23_sub'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U2'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];

dir = '';

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_SUB_refCA1_forAnalysis' '.xlsx']);
ReactTable.SUB = readtable([ROOT.Save '\ReactTable_SUB_SUB.xlsx']);
UnitsTable.SUB = readtable([ROOT.Save '\UnitsTable_SUB_forAnalysis.xlsx']);
UnitsTable_field.SUB = readtable([ROOT.Save '\UnitsTable_SUB_field_forAnalysis.xlsx']);

RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx']);
ReactTable.CA1 = readtable([ROOT.Save '\ReactTable_CA1_CA1.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_forAnalysis.xlsx']);
UnitsTable_field.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx']);

%%

thisRegion = 'SUB'; cl=0;

UT = UnitsTable.(thisRegion);
UTf = UnitsTable_field.(thisRegion); 
RT = RipplesTable.(thisRegion);
ReT = ReactTable.(thisRegion);

UnitPair=table; u=1;
for uid=1:size(UT,1)
    p0 = sum(RT.rat==UT.rat(uid) & RT.session==UT.session(uid));
    thisUnits = UT(UT.rat==UT.rat(uid) & UT.session==UT.session(uid),:);
    thisRips = unique(ReT(find(cellfun(Params.cellfind(UT.ID{uid}),ReT.UnitID)),[1 2]));
    p1 = size(thisRips,1);
    for uid2=1:size(thisUnits,1)
        if strncmp(thisUnits.ID{uid2},UT.ID{uid},12), continue; end

        thisRips2=unique(ReT(find(cellfun(Params.cellfind(thisUnits.ID{uid2}),ReT.UnitID)),[1 2]));
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

        thisField_1 = UTf(find(strncmp(UTf.ID,UT.ID{uid},12)),:);
        thisField_2 = UTf(find(strncmp(UTf.ID,UT.ID{uid2},12)),:);
        if isempty(thisField_1) | isempty(thisField_2)
            UnitPair.region(u)=1;
            continue;
        end

        UnitPair.L1(u)={thisField_1.RDI_LScene};
        UnitPair.R1(u)={thisField_1.RDI_RScene};
        UnitPair.C1(u)={thisField_1.RDI_LR};

        UnitPair.L2(u)={thisField_2.RDI_LScene};
        UnitPair.R2(u)={thisField_2.RDI_RScene};
        UnitPair.C2(u)={thisField_2.RDI_LR};
%%
        UT.Selectivity_LScene(uid) = getSel(UnitPair.L1{u});
        UT.Selectivity_RScene(uid) = getSel(UnitPair.R1{u});
         UT.Selectivity_LR(uid) = getSel(UnitPair.C1{u});

                 thisUnits.Selectivity_LScene(uid2) = getSel(UnitPair.L2{u});
        thisUnits.Selectivity_RScene(uid2) = getSel(UnitPair.R2{u});
         thisUnits.Selectivity_LR(uid2) = getSel(UnitPair.C2{u});
%%
        UnitPair.Lt1(u) = UT.Selectivity_LScene(uid) * max(UnitPair.L1{u}) / max([10^(-300),abs(max(UnitPair.L1{u}))]);
        UnitPair.Rt1(u) = UT.Selectivity_RScene(uid) * max(UnitPair.R1{u}) / max([10^(-300),abs(max(UnitPair.R1{u}))]);
        UnitPair.Ct1(u) = UT.Selectivity_LR(uid) * max(UnitPair.C1{u}) / max([10^(-300),abs(max(UnitPair.C1{u}))]);

        UnitPair.Lt2(u) = thisUnits.Selectivity_LScene(uid2) * max(UnitPair.L2{u}) / max([10^(-300),abs(max(UnitPair.L2{u}))]);
        UnitPair.Rt2(u) = thisUnits.Selectivity_RScene(uid2) * max(UnitPair.R2{u}) / max([10^(-300),abs(max(UnitPair.R2{u}))]);
        UnitPair.Ct2(u) = thisUnits.Selectivity_LR(uid2) * max(UnitPair.C2{u}) / max([10^(-300),abs(max(UnitPair.C2{u}))]);

        UnitPair.Lm1(u)=mean(thisField_1.RDI_LScene(abs(thisField_1.RDI_LScene)>0.1));
        UnitPair.Rm1(u)=mean(thisField_1.RDI_RScene(abs(thisField_1.RDI_RScene)>0.1));
        UnitPair.Cm1(u)=mean(thisField_1.RDI_LR(abs(thisField_1.RDI_LR)>0.1));

      UnitPair.Lm2(u)=mean(thisField_2.RDI_LScene(abs(thisField_2.RDI_LScene)>0.1));
        UnitPair.Rm2(u)=mean(thisField_2.RDI_RScene(abs(thisField_2.RDI_RScene)>0.1));
        UnitPair.Cm2(u)=mean(thisField_2.RDI_LR(abs(thisField_2.RDI_LR)>0.1));


        UnitPair.region(u) = cl;
%              UnitPair.NumField1(u) = UT.NumField(uid);
%              UnitPair.NumField2(u) = thisUnits.NumField(uid2);

%         UnitPair.NumField1(u) = sum(UT.AvgFR(uid)==UT.AvgFR);
%          UnitPair.NumField2(u) = sum(thisUnits.AvgFR(uid2)==thisUnits.AvgFR);
        
        u=u+1;
    end
end
writetable(UnitPair,[ROOT.Save '\UnitPair_' thisRegion '.xlsx'],'writemode','replacefile');

%%
UnitPair_CA1 = readtable([ROOT.Save '\UnitPair_CA1.xlsx']);
UnitPair_SUB = readtable([ROOT.Save '\UnitPair_SUB.xlsx']);
%%
UnitPair = UnitPair_SUB;
UnitPair.CoactP = UnitPair.pq./UnitPair.un;
CR_SUB.SF = UnitPair.CoactP((abs(UnitPair.Lt1)==1 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==1 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==1 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));
CR_SUB.SMF = UnitPair.CoactP((abs(UnitPair.Lt1)==2 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==2 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==2 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));
CR_SUB.MF = UnitPair.CoactP((abs(UnitPair.Lt1)==3 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==3 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==3 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));

UnitPair = UnitPair_CA1;
UnitPair.CoactP = UnitPair.pq./UnitPair.un;
CR_CA1.SF = UnitPair.CoactP((abs(UnitPair.Lt1)==1 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==1 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==1 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));
CR_CA1.SMF = UnitPair.CoactP((abs(UnitPair.Lt1)==2 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==2 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==2 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));
CR_CA1.MF = UnitPair.CoactP((abs(UnitPair.Lt1)==3 & abs(UnitPair.Lt1)== abs(UnitPair.Lt2)) |...
    (abs(UnitPair.Rt1)==3 & abs(UnitPair.Rt1)== abs(UnitPair.Rt2)) |...
    (abs(UnitPair.Ct1)==3 & abs(UnitPair.Ct1)== abs(UnitPair.Ct2)));


data = [nanmean([CR_SUB.SF;CR_CA1.SF]) nanmean(CR_SUB.SF) nanmean(CR_CA1.SF);...
    nanmean([CR_SUB.SMF;CR_SUB.MF;CR_CA1.SMF;CR_CA1.MF]) nanmean([CR_SUB.SMF;CR_SUB.MF]) nanmean([CR_CA1.SMF;CR_CA1.MF])];

bar(data')
xticklabels({'All','SUB','CA1'})
legend({'SF-SF','MF-MF'},'location','best')
ylabel('Co-Reactivation Probability')
set(gca,'fontsize',12,'fontweight','b')
[h,p] = ttest2([CR_SUB.SF;CR_CA1.SF],[CR_SUB.SMF;CR_SUB.MF;CR_CA1.SMF;CR_CA1.MF])
%%
UnitPair = UnitPair_SUB;
UnitPair.CoactP = UnitPair.pq./UnitPair.un;
d = 'C';
Coreacts.([d '4']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])>0,:);
Coreacts.([d '4']) = Coreacts.([d '4']).CoactP;
Coreacts.([d '3']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])<0,:);
Coreacts.([d '3']) = Coreacts.([d '3']).CoactP;
Coreacts.([d '2']) = [UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))>=1,:);...
    UnitPair(abs(UnitPair.([d 't2']))<1 & abs(UnitPair.([d 't1']))>=1,:)];
Coreacts.([d '2']) = Coreacts.([d '2']).CoactP;
Coreacts.([d '1']) = UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))<1,:);
Coreacts.([d '1']) = Coreacts.([d '1']).CoactP;

data_sub = [nanmean(Coreacts.L1) nanmean(Coreacts.R1) nanmean(Coreacts.C1);...
    nanmean(Coreacts.L2) nanmean(Coreacts.R2) nanmean(Coreacts.C2);...
    nanmean(Coreacts.L3) nanmean(Coreacts.R3) nanmean(Coreacts.C3);...
    nanmean(Coreacts.L4) nanmean(Coreacts.R4) nanmean(Coreacts.C4)];
subplot(1,2,1)
bar(data_sub')
ylabel('Co-Reactivation Probability'); ylim([0 0.16])
set(gca,'fontsize',12,'fontweight','b')

UnitPair = UnitPair_CA1;
UnitPair.CoactP = UnitPair.pq./UnitPair.un;
d = 'C';
Coreacts.([d '4']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])>0,:);
Coreacts.([d '4']) = Coreacts.([d '4']).CoactP;
Coreacts.([d '3']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])<0,:);
Coreacts.([d '3']) = Coreacts.([d '3']).CoactP;
Coreacts.([d '2']) = [UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))>=1,:);...
    UnitPair(abs(UnitPair.([d 't2']))<1 & abs(UnitPair.([d 't1']))>=1,:)];
Coreacts.([d '2']) = Coreacts.([d '2']).CoactP;
Coreacts.([d '1']) = UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))<1,:);
Coreacts.([d '1']) = Coreacts.([d '1']).CoactP;

data_ca1 = [nanmean(Coreacts.L1) nanmean(Coreacts.R1) nanmean(Coreacts.C1);...
    nanmean(Coreacts.L2) nanmean(Coreacts.R2) nanmean(Coreacts.C2);...
    nanmean(Coreacts.L3) nanmean(Coreacts.R3) nanmean(Coreacts.C3);...
    nanmean(Coreacts.L4) nanmean(Coreacts.R4) nanmean(Coreacts.C4)];
subplot(1,2,2)
bar(data_ca1')
ylabel('Co-Reactivation Probability'); ylim([0 0.16])
set(gca,'fontsize',12,'fontweight','b')
% data_ca1 = [nanmean([Coreacts.L1; Coreacts.R1;Coreacts.C1]);...
%     nanmean([Coreacts.L2; Coreacts.R2;Coreacts.C2]);...
%     nanmean([Coreacts.L3; Coreacts.R3;Coreacts.C3]);...
%     nanmean([Coreacts.L4; Coreacts.R4;Coreacts.C4])];

% [h,p] = ttest2(Coreacts.L4,Coreacts.L4)
%%
subplot(1,2,1)
UP = UnitPair_SUB(abs(UnitPair_SUB.Lt1)>=1 & abs(UnitPair_SUB.Lt2)>=1,:);
UP.CoactP = UP.pq./UP.un;
x = abs(UP.Lm1-UP.Lm2); y=UP.CoactP; 
y(isnan(x))=[]; x(isnan(x))=[];

scatter(x,y,20,'k')
xlabel('Left scene selectivity difference'); xlim([0 1.8])
ylabel('Co-Reactivation Probability'); ylim([0 1])

set(gca,'fontsize',12,'fontweight','b')

[r,p] = corr(x,y)

subplot(1,2,2)
UP = UnitPair_CA1(abs(UnitPair_CA1.Lt1)>=1 & abs(UnitPair_CA1.Lt2)>=1,:);
UP.CoactP = UP.pq./UP.un;
x = abs(UP.Lm1-UP.Lm2); y=UP.CoactP; 
y(isnan(x))=[]; x(isnan(x))=[];
x(isnan(y))=[]; y(isnan(y))=[];

scatter(x,y,20,'k')
xlabel('Left scene selectivity difference'); xlim([0 1.8])
ylabel('Co-Reactivation Probability'); ylim([0 1])
set(gca,'fontsize',12,'fontweight','b')
[r,p] = corr(x,y)


%%
function s = getSel(rdis)
rdisig = rdis(abs(rdis)>0.1);
if length(rdis)>1
    if length(rdisig)==1
        s=2;
    elseif length(rdisig)>1
        if min(rdisig>0) | min(rdisig<0)
            s=3;
        else
            s=4;
        end
    else
        s=0.5;
    end

else
    if length(rdisig)>0
        s=1;
    else
        s=0;
    end

end
end



