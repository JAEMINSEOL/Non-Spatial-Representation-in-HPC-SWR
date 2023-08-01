Initial_SWRFilter_common;
warning off

ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip5 = [ROOT.Mother '\Processed Data\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R23_sub'];
ROOT.Unit1 = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U2'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];

dir = '';

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
RegionList = {'CA1','SUB'};

RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_SUB_refCA1_forAnalysis' '.xlsx']);
ReactTable.SUB = readtable([ROOT.Save '\ReactTable_SUB_SUB.xlsx']);
UnitsTable.SUB = readtable([ROOT.Save '\UnitsTable_SUB_forAnalysis_TP.xlsx']);
UnitsTable_field.SUB = readtable([ROOT.Save '\UnitsTable_SUB_field_forAnalysis.xlsx']);

RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx']);
ReactTable.CA1 = readtable([ROOT.Save '\ReactTable_CA1_CA1.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_forAnalysis_TP.xlsx']);
UnitsTable_field.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx']);

FRMaps.SUB = LoadFRMap(ROOT,UnitsTable.SUB);
FRMaps.SUB_field = LoadFRMap(ROOT,UnitsTable_field.SUB);
FRMaps.CA1 = LoadFRMap(ROOT,UnitsTable.CA1);
FRMaps.CA1_field = LoadFRMap(ROOT,UnitsTable_field.CA1);
%%
for Reg = 1:numel(RegionList)
    thisRegion = RegionList{Reg};
    % thisRegion = 'SUB';
    cl=Reg;

    UT = UnitsTable.(thisRegion);
    UTf = UnitsTable_field.(thisRegion);
    RT = RipplesTable.(thisRegion);
    ReT = ReactTable.(thisRegion);
    FR  = FRMaps.(thisRegion);

    UnitPair=table; u=1;
    for uid=1:size(UT,1)
        p0 = sum(RT.rat==UT.rat(uid) & RT.session==UT.session(uid));
        thisUnits = UT(UT.rat==UT.rat(uid) & UT.session==UT.session(uid),:);
        thisRips = unique(ReT(find(cellfun(Params.cellfind(UT.ID{uid}),ReT.UnitID)),[1 2]));
        p1 = size(thisRips,1);

        FR_u1 = FR(1,:,uid);
        FR_u1(isnan(FR_u1)) = 0;
        FR_u1 = FR_u1 ./ max(FR_u1);

        for uid2=1:size(thisUnits,1)
            if strncmp(thisUnits.ID{uid2},UT.ID{uid},12), continue; end

            if uid~=1
                if(find(strcmp(UT.ID{uid},UnitPair.UID2) & strcmp(thisUnits.ID{uid2},UnitPair.UID1)))
                    continue;
                end
            end
            thisRips2=unique(ReT(find(cellfun(Params.cellfind(thisUnits.ID{uid2}),ReT.UnitID)),[1 2]));
            CoActRips = intersect(thisRips.RippleID,thisRips2.RippleID);

            uid1_2 = find(strcmp(thisUnits.ID{uid2},UT.ID));

            FR_u2 = FR(1,:,uid1_2);
            FR_u2(isnan(FR_u2)) = 0;
            FR_u2 = FR_u2 ./ max(FR_u2);

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
            thisField_2 = UTf(find(strncmp(UTf.ID,thisUnits.ID{uid2},12)),:);
            if isempty(thisField_1) | isempty(thisField_2)
                UnitPair.region(u)=1;
                continue;
            end

            RDI_L1 = []; RDI_R1 = []; RDI_C1 = [];
            for f = 1:size(thisField_1,1)
                F = load([ROOT.Unit1 '\' thisField_1.ID{f} '.mat'],'RDIs_field'); F= F.RDIs_field;
                RDI_L1 = [RDI_L1;interp1(1:size(F,2),F(1,:),linspace(1,size(F,2),40))];
                RDI_R1 = [RDI_R1; interp1(1:size(F,2),F(2,:),linspace(1,size(F,2),40))];
                RDI_C1 = [RDI_C1; interp1(1:size(F,2),F(3,:),linspace(1,size(F,2),40))];
            end

            [~,o] = sort(thisField_1.RDI_LScene); RDI_L1 = RDI_L1(o,:);
            [~,o] = sort(thisField_1.RDI_RScene); RDI_R1 = RDI_R1(o,:);
            [~,o] = sort(thisField_1.RDI_LR); RDI_C1 = RDI_C1(o,:);


             RDI_L2 = []; RDI_R2 = []; RDI_C2 = [];
            for f = 1:size(thisField_2,1)
                F = load([ROOT.Unit1 '\' thisField_2.ID{f} '.mat'],'RDIs_field'); F= F.RDIs_field;
                RDI_L2 = [RDI_L2;interp1(1:size(F,2),F(1,:),linspace(1,size(F,2),40))];
                RDI_R2 = [RDI_R2; interp1(1:size(F,2),F(2,:),linspace(1,size(F,2),40))];
                RDI_C2 = [RDI_C2; interp1(1:size(F,2),F(3,:),linspace(1,size(F,2),40))];
            end

            [~,o] = sort(thisField_2.RDI_LScene); RDI_L2 = RDI_L2(o,:);
            [~,o] = sort(thisField_2.RDI_RScene); RDI_R2 = RDI_R2(o,:);
            [~,o] = sort(thisField_2.RDI_LR); RDI_C2 = RDI_C2(o,:);

            s1 = size(RDI_L1,1); s2 = size(RDI_L2,1);
            if s1<s2,RDI_L1(s1+1:s2,:)=nan; elseif s1>s2,RDI_L2(s2+1:s1,:)=nan; end
            RDI_L1(isnan(RDI_L1))=0; RDI_L2(isnan(RDI_L2))=0;
           

             s1 = size(RDI_R1,1); s2 = size(RDI_R2,1);
            if s1<s2,RDI_R1(s1+1:s2,:)=nan; elseif s1>s2,RDI_R2(s2+1:s1,:)=nan; end
            RDI_R1(isnan(RDI_R1))=0; RDI_R2(isnan(RDI_R2))=0;
            

             s1 = size(RDI_C1,1); s2 = size(RDI_C2,1);
            if s1<s2,RDI_C1(s1+1:s2,:)=nan; elseif s1>s2,RDI_C2(s2+1:s1,:)=nan; end
            RDI_C1(isnan(RDI_C1))=0; RDI_C2(isnan(RDI_C2))=0;


            UnitPair.Sp(u) = corr(FR_u1',FR_u2');

            UnitPair. Nsp_L(u) = corr2(RDI_L1,RDI_L2);
            UnitPair. Nsp_R(u) = corr2(RDI_R1,RDI_R2);
            UnitPair. Nsp_C(u) = corr2(RDI_C1,RDI_C2);

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

            UnitPair.Lm1(u)=mean(thisField_1.RDI_LScene(abs(thisField_1.RDI_LScene)>0));
            UnitPair.Rm1(u)=mean(thisField_1.RDI_RScene(abs(thisField_1.RDI_RScene)>0));
            UnitPair.Cm1(u)=mean(thisField_1.RDI_LR(abs(thisField_1.RDI_LR)>0));

            UnitPair.Lm2(u)=mean(thisField_2.RDI_LScene(abs(thisField_2.RDI_LScene)>0));
            UnitPair.Rm2(u)=mean(thisField_2.RDI_RScene(abs(thisField_2.RDI_RScene)>0));
            UnitPair.Cm2(u)=mean(thisField_2.RDI_LR(abs(thisField_2.RDI_LR)>0));


            UnitPair.region(u) = cl;
            %              UnitPair.NumField1(u) = UT.NumField(uid);
            %              UnitPair.NumField2(u) = thisUnits.NumField(uid2);

            %         UnitPair.NumField1(u) = sum(UT.AvgFR(uid)==UT.AvgFR);
            %          UnitPair.NumField2(u) = sum(thisUnits.AvgFR(uid2)==thisUnits.AvgFR);

            if u==1032
                1;
            end

            u=u+1;
        end
    end
     save([ROOT.Save '\UnitPair_' thisRegion '.mat'],'UnitPair');
end
% writetable(UnitPair,[ROOT.Save '\UnitPair_' thisRegion '.xlsx'],'writemode','replacefile');

%%
UnitPairT=struct;
UnitPairT.CA1 = load([ROOT.Save '\UnitPair_CA1.mat']); UnitPairT.CA1= UnitPairT.CA1.UnitPair;
UnitPairT.SUB = load([ROOT.Save '\UnitPair_SUB.mat']); UnitPairT.SUB= UnitPairT.SUB.UnitPair;
%%
UnitPair = UnitPairT.SUB;
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

UnitPair = UnitPairT.CA1;
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

figure;
bar(data')
xticklabels({'All','SUB','CA1'})
legend({'SF-SF','MF-MF'},'location','best')
ylabel('Co-Reactivation Probability')
set(gca,'fontsize',12,'fontweight','b')
[h,p] = ttest2([CR_CA1.SF],[CR_CA1.SMF;CR_CA1.MF])
%%  SF-SF, MF-MF coreactivation prob.
figure;
UnitPair = UnitPairT.SUB;
% UnitPair(UnitPair.Lt1==1 | UnitPair.Lt1==0 | UnitPair.Lt2==1 | UnitPair.Lt2==0,:) = [];
UnitPair(UnitPair.Lt1>1 | UnitPair.Lt1==0.5 | UnitPair.Lt2>1 | UnitPair.Lt2==0.5,:) = [];
UnitPair.CoactP = UnitPair.pq./UnitPair.un;
d_List = {'L','R','C'};

for d0 = 1:numel(d_List)
d = d_List{d0};
Coreacts.([d '4']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])>0,:);
Coreacts.([d '4']) = Coreacts.([d '4']).CoactP;
Coreacts.([d '3']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])<0,:);
Coreacts.([d '3']) = Coreacts.([d '3']).CoactP;
Coreacts.([d '2']) = [UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))>=1,:);...
    UnitPair(abs(UnitPair.([d 't2']))<1 & abs(UnitPair.([d 't1']))>=1,:)];
Coreacts.([d '2']) = Coreacts.([d '2']).CoactP;
Coreacts.([d '1']) = UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))<1,:);
Coreacts.([d '1']) = Coreacts.([d '1']).CoactP;
end

data_sub = [nanmean(Coreacts.L1) nanmean(Coreacts.R1) nanmean(Coreacts.C1);...
    nanmean(Coreacts.L2) nanmean(Coreacts.R2) nanmean(Coreacts.C2);...
    nanmean(Coreacts.L3) nanmean(Coreacts.R3) nanmean(Coreacts.C3);...
    nanmean(Coreacts.L4) nanmean(Coreacts.R4) nanmean(Coreacts.C4)];
subplot(1,4,2)
bar(data_sub')
ylabel('Co-Reactivation Probability'); ylim([0 0.18])
set(gca,'fontsize',12,'fontweight','b')
title('SUB')
xticklabels({'Left','Right','Choice'})


UnitPair = UnitPairT.CA1;
% UnitPair(UnitPair.Lt1==1 | UnitPair.Lt1==0 | UnitPair.Lt2==1 | UnitPair.Lt2==0,:) = [];
UnitPair(UnitPair.Lt1>1 | UnitPair.Lt1==0.5 | UnitPair.Lt2>1 | UnitPair.Lt2==0.5,:) = [];
UnitPair.CoactP = UnitPair.pq./UnitPair.un;

for d0 = 1:numel(d_List)
d = d_List{d0};
Coreacts.([d '4']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])>0,:);
Coreacts.([d '4']) = Coreacts.([d '4']).CoactP;
Coreacts.([d '3']) = UnitPair(abs(UnitPair.([d 't1']))>=1 & abs(UnitPair.([d 't2']))>=1 & UnitPair.([d 't1']).*UnitPair.([d 't2'])<0,:);
Coreacts.([d '3']) = Coreacts.([d '3']).CoactP;
Coreacts.([d '2']) = [UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))>=1,:);...
    UnitPair(abs(UnitPair.([d 't2']))<1 & abs(UnitPair.([d 't1']))>=1,:)];
Coreacts.([d '2']) = Coreacts.([d '2']).CoactP;
Coreacts.([d '1']) = UnitPair(abs(UnitPair.([d 't1']))<1 & abs(UnitPair.([d 't2']))<1,:);
Coreacts.([d '1']) = Coreacts.([d '1']).CoactP;
end

data_ca1 = [nanmean(Coreacts.L1) nanmean(Coreacts.R1) nanmean(Coreacts.C1);...
    nanmean(Coreacts.L2) nanmean(Coreacts.R2) nanmean(Coreacts.C2);...
    nanmean(Coreacts.L3) nanmean(Coreacts.R3) nanmean(Coreacts.C3);...
    nanmean(Coreacts.L4) nanmean(Coreacts.R4) nanmean(Coreacts.C4)];
subplot(1,2,2)
bar(data_ca1')
ylabel('Co-Reactivation Probability'); ylim([0 0.18])
set(gca,'fontsize',12,'fontweight','b')
title('CA1')
xticklabels({'Left','Right','Choice'})
legend({'Non-Non','Non-Remap','Remap-Remap (diff.)','Remap-Remap (same)'},'location','eastoutside','Orientation','vertical')
% data_ca1 = [nanmean([Coreacts.L1; Coreacts.R1;Coreacts.C1]);...
%     nanmean([Coreacts.L2; Coreacts.R2;Coreacts.C2]);...
%     nanmean([Coreacts.L3; Coreacts.R3;Coreacts.C3]);...
%     nanmean([Coreacts.L4; Coreacts.R4;Coreacts.C4])];

 [h,p] = ttest2(Coreacts.C4,Coreacts.C3)

 %%
 data = [data_sub(:,1)' ; data_ca1(:,1)'];
figure;
bar(data)
ylabel('Co-Reactivation Probability');
set(gca,'fontsize',12,'fontweight','b')
title('SF')
xticklabels({'SUB','CA1'})
legend({'Non-Non','Non-Remap','Remap-Remap (diff.)','Remap-Remap (same)'},'location','eastoutside','Orientation','vertical')

%%
figure;
d='L';
for r = 1:numel(RegionList)
    thisRegion = RegionList{r};
figure; sgtitle(thisRegion)
for i=1:3
UP = UnitPairT.(thisRegion);
UP.CoactP = UP.pq./UP.un;
UP = UP(abs(UP.Sp)>1-i*0.3 & abs(UP.Sp)<=1.3-i*0.3,:);
x = 1-abs(UP.([d 'm1'])-UP.([d 'm2'])); s=UP.Sp; y=UP.CoactP;
id1 = ((UP.([d 't1'])>1|UP.([d 't1'])==0.5)& UP.([d 't1'])<4) & ((UP.([d 't2'])>1|UP.([d 't2'])==0.5)& UP.([d 't1'])<4);
id2 = (UP.([d 't1'])==1|UP.([d 't1'])==0) & (UP.([d 't2'])==1|UP.([d 't2'])==0);

ax = subplot(2,3,i);
% y(isnan(x))=[]; x(isnan(x))=[];
scatter(x(id1),y(id1),40,'k')
% scatter3(x(id2),s(id2),y(id2),20,'b')
xlabel('Non-Spatial Similiarity'); xlim([0 1])
ylabel('Co-Reactivation Probability'); ylim([0 1])
% zlabel('Co-Reactivation Probability'); zlim([0 1])
title(['MF, ' num2str(1-i*0.3) ' < Sp.Sim. < ' num2str(1.3-i*0.3)])
set(gca,'fontsize',12,'fontweight','b')
lsline(ax)
text(0.5,.8,['corr = ' jjnum2str(corr(x(id1),y(id1)),3)])

ax = subplot(2,3,i+3);
% y(isnan(x))=[]; x(isnan(x))=[];
scatter(x(id2),y(id2),40,'k')
% scatter3(x(id2),s(id2),y(id2),20,'b')
xlabel('Non-Spatial Similiarity'); xlim([0 1])
ylabel('Co-Reactivation Probability'); ylim([0 1])
% zlabel('Co-Reactivation Probability'); zlim([0 1])
title(['SF, ' num2str(1-i*0.3) ' < Sp.Sim. < ' num2str(1.3-i*0.3)])
set(gca,'fontsize',12,'fontweight','b')
lsline(ax)

x1 = x(id2); y1 = y(id2);
y1(isnan(x1))=[]; x1(isnan(x1))=[]; 
x1(isnan(y1))=[]; y1(isnan(y1))=[]; 
text(0.5,.8,['corr = ' jjnum2str(corr(x1,y1),3)])


end
end


%%
UP = UnitPairT.SUB;
UP = UP((UP.Lt1>1 | UP.Lt1==0.5) & (UP.Lt2>1 | UP.Lt2==0.5),:);

% UP = UP((UP.Lt1==1 | UP.Lt1==0) & (UP.Lt2==1 | UP.Lt2==0),:);
UP.CoactP = UP.pq ./ UP.un;
x2 = 1-abs(UP.([d 'm1'])-UP.([d 'm2'])); x1=abs(UP.Sp); y=UP.CoactP;
X = [ones(size(x1)) x1 x2 x1.*x2];
b = regress(y,X)    % Removes NaN data

scatter3(x1,x2,y,'filled')
hold on
x1fit = min(x1):.01:max(x1);
x2fit = min(x2):.01:max(x2);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT + b(4)*X1FIT.*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('Sp Sim.')
ylabel('NonSp Sim.')
zlabel('CoReact Prob')
view(50,10)
hold off
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



