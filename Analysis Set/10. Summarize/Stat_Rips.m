Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Save '\ripples_mat\R4'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R11_sub'];
ROOT.Units = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];



thisRegion0 = 'CA1';
thisRegion = 'CA1';
thisRegion2 = [thisRegion '_field'];
RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis' '.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable.CA1_field  = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);

thisRegion0 = 'SUB';
thisRegion = 'SUB';
thisRegion2 = [thisRegion '_field'];
RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis' '.xlsx']);
UnitsTable.SUB = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable.SUB_field = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);

CList = [ [207 8 23]/255;[23 84 181]/255];

RipplesTable.CA1.NS = nanmin([RipplesTable.CA1.pRDI_L_UV,RipplesTable.CA1.pRDI_R_UV,RipplesTable.CA1.pRDI_C_UV],[],2) <0.05;
RipplesTable.CA1.NS_p = nanmin([RipplesTable.CA1.pRDI_L_UV,RipplesTable.CA1.pRDI_R_UV,RipplesTable.CA1.pRDI_C_UV],[],2);

RipplesTable.SUB.NS = nanmin([RipplesTable.SUB.pRDI_L_UV,RipplesTable.SUB.pRDI_R_UV,RipplesTable.SUB.pRDI_C_UV],[],2) <0.05;
RipplesTable.SUB.NS_p = nanmin([RipplesTable.SUB.pRDI_L_UV,RipplesTable.SUB.pRDI_R_UV,RipplesTable.SUB.pRDI_C_UV],[],2);

RipplesTable.CA1.S = RipplesTable.CA1.DecodingP_all <0.05;
RipplesTable.CA1.S_p = RipplesTable.CA1.DecodingP_all ;

RipplesTable.SUB.S = RipplesTable.SUB.DecodingP_all  <0.05;
RipplesTable.SUB.S_p = RipplesTable.SUB.DecodingP_all ;

RipplesTable.CA1.NS(RipplesTable.CA1.nFields<5)=0;
RipplesTable.CA1.NS_p(RipplesTable.CA1.nFields<5)=nan;

RipplesTable.SUB.NS(RipplesTable.SUB.nFields<5)=0;
RipplesTable.SUB.NS_p(RipplesTable.SUB.nFields<5)=nan;

T0 = RipplesTable.SUB;
T1 = RipplesTable.CA1;

%% spike, ensemble 다시 계산 
thisRegion = 'SUB'; cl=1;
thisRegion2 = [thisRegion '_field'];
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);

for rid = 1:size(T0,1)
    thisReact = ReactTable(find(cellfun(Params.cellfind(T0.ID{rid}),ReactTable.RippleID)),:);

    RipplesTable.(thisRegion).spike(rid) = size(thisReact,1);
    RipplesTable.(thisRegion).ensemble(rid) = size(unique(thisReact.UnitID),1);
end

% writetable(RipplesTable.(thisRegion), [ROOT.Save '\RipplesTable_' thisRegion2 '_RDIs_UV_cell_HeteroIn_AllPopul' '.xlsx'],'WriteMode','replacefile');



%% SUB vs. CA1 spike, ensemble 비교

x1 = T0.spike ./ T0.ensemble;
x2 = T1.spike ./ T1.ensemble;
leg = {'SUB','CA1'};
HistAndBar(x1,x2,CList,leg,'spikes per cell')


[h,p,c,stats] = ttest2(x1,x2)
%%

figure; hold on
T0 = RipplesTable.SUB;
T1 = RipplesTable.CA1;


histogram(T0.NS_p,'binwidth',0.01,'Normalization','probability','FaceColor',CList(1,:))
histogram(T1.NS_p,'binwidth',0.01,'Normalization','probability','FaceColor',CList(2,:))

legend({'SUB','CA1'})
set(gca,'fontsize',12,'fontweight','b')
xlabel('p-value for non-spatial selectivity permutation test')
ylabel('proportion')
%%
CList = [[139 193 69]/255; [29 111 169]/255];
figure; hold on
T0 = RipplesTable.SUB;
T1 = RipplesTable.CA1;
histogram(T0.S_p,'binwidth',0.01,'Normalization','probability','FaceColor',CList(1,:))
histogram(T1.S_p,'binwidth',0.01,'Normalization','probability','FaceColor',CList(2,:))

legend({'SUB','CA1'})
set(gca,'fontsize',12,'fontweight','b')
xlabel('p-value for position decoding permutation test')
ylabel('proportion')
%%
figure; hold on
x = [1 2];
data = [nanmean(T0.NS_p) nanmean(T1.NS_p); nanmean(T0.S_p) nanmean(T1.S_p)];
err = [nanstd(T0.NS_p)/sqrt(size(T0,1)) nanstd(T1.NS_p)/sqrt(size(T1,1));...
    nanstd(T0.S_p)/sqrt(size(T0,1)) nanstd(T1.S_p)/sqrt(size(T1,1))];

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
%% 
figure; hold on
T0 = RipplesTable.SUB;
T1 = RipplesTable.CA1;
c0= cdfplot(T0.NS_p); c0.Color = CList(1,:);
c1 = cdfplot(T1.NS_p); c1.Color = CList(2,:);

%%
sum(T0.NS_p<0.05 & ~(T0.S_p<0.05))/size(T0,1)
sum(~(T0.NS_p<0.05) & (T0.S_p<0.05))/size(T0,1)

sum(T1.NS_p<0.05 & ~(T1.S_p<0.05))/size(T1,1)
sum(~(T1.NS_p<0.05) & (T1.S_p<0.05))/size(T1,1)

sum((T0.NS_p<0.05) & (T0.S_p<0.05))/size(T0,1)
sum(T1.NS_p<0.05 & (T1.S_p<0.05))/size(T1,1)
%%
[p,h] = ttest2(T0.NS_p,T1.NS_p)

%% SF vs. MF 
thisRegion = 'SUB'; cl=1;
RT = RipplesTable.(thisRegion);
UnitRipTable_Field = [];
thisRegion2 = [thisRegion '_field'];
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);
UT = UnitsTable.(thisRegion2); UT = UT(UT.NumField>0,:);

for rid=1:size(RT,1)
     thisRip =RT(rid,:);
         thisReact_A = ReactTable(strcmp(ReactTable.RippleID,thisRip.ID{1}),:);
           thisUnits=table;
    for u=1:size(thisReact_A,1)
        thisUnits =[thisUnits;UT(find(cellfun(Params.cellfind(thisReact_A.UnitID(u)),UT.ID)),:)];
    end

    UnitRipTable_Field(rid,1) = sum(thisUnits.NumField==1)/size(thisUnits,1);
    UnitRipTable_Field(rid,2) = sum(thisUnits.NumField==2)/size(thisUnits,1);
    UnitRipTable_Field(rid,3) = sum(thisUnits.NumField==3)/size(thisUnits,1);
    UnitRipTable_Field(rid,4) = sum(thisUnits.NumField>3)/size(thisUnits,1);
end

%% plot
figure; hold on

x = [1 2];
x1 = UnitRipTable_Field(:,1);
x2 = sum(UnitRipTable_Field(:,2:4),2);

data = [nanmean(x1) nanmean(x2)];
err = [nanstd(x1)/sqrt(size(x1,1)) nanstd(x2)/sqrt(size(x2,1))];


b = bar(x,data);
b.FaceColor='flat';
b.CData(1,:) = [1 1 1];
b.CData(2,:) = CList(cl,:);
b.FaceAlpha = 0.3;

    % get x positions per group
    xpos = b(1).XData + b(1).XOffset;
    % draw errorbar
    errorbar(xpos, data(1,:), err(1,:), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);

    t1 = sum(UT.NumField==1)/size(UT,1);
    t2 = sum(UT.NumField>1)/size(UT,1);

%      hold on
%     c = bar(x,[t1 t2]); 
%     c.FaceAlpha=0;
%     c.EdgeColor='r';
% 
%     [~,p1] = ttest(x1,t1);
%     [~,p2] = ttest(x2,t2);
% 
%     if p1>0.0001, text(1-.2,data(1)+.05,jjnum2str(p1,4)); else, text(1-.2,data(1)+.05,'p<0.0001'); end
%     if p2>0.0001, text(2-.2,data(2)+.05,jjnum2str(p2,4)); else, text(2-.2,data(2)+.05,'p<0.0001'); end

ylim([0 1])
    xticks([1 2]); xticklabels({'SF','MF'})
%% Unit type Rip table
UnitRipTable = table('Size',[6 17], 'VariableTypes', {'string' 'double' 'double' 'double' 'double',...
    'double' 'double' 'double' 'double' 'double' 'double' 'double' 'double' 'double' 'double' 'double' 'double'},...
    'VariableNames',{'Unit','X_X','Sp_X',...
    'Sp_NSp_L' 'Sp_NSp_R' 'Sp_NSp_C' 'Sp_NSp_L_R' 'Sp_NSp_L_C' 'Sp_NSp_R_C' 'Sp_NSp_L_R_C',...
    'X_NSp_L' 'X_NSp_R' 'X_NSp_C' 'X_NSp_L_R' 'X_NSp_L_C' 'X_NSp_R_C' 'X_NSp_L_R_C'});
UnitRipTable.Unit = {'SF_N'; 'SF_S' ; 'MF_N' ; 'MF_S' ; 'MF_Homo'; 'MF_Hetero'};

UnitRipTable_L = UnitRipTable; UnitRipTable_R = UnitRipTable; UnitRipTable_C = UnitRipTable;

for rid=1:size(RipplesTable,1)
     thisRip =RipplesTable(rid,:);
         thisReact_A = ReactTable(strcmp(ReactTable.RippleID,thisRip.ID{1}),:);
           thisUnits=table;
    for u=1:size(thisReact_A,1)
        thisUnits =[thisUnits;UnitsTable_A(find(cellfun(Params.cellfind(thisReact_A.UnitID(u)),UnitsTable_A.ID)),:)];
    end
    rsuf = thisRip.Type{1};

    for uid = 1:size(thisUnits,1)

        switch thisUnits.Selectivity_LScene(uid)
            case 0, UnitRipTable_L.(rsuf)(1) = UnitRipTable_L.(rsuf)(1)+1;
            case 1, UnitRipTable_L.(rsuf)(2) = UnitRipTable_L.(rsuf)(2)+1;
            case 0.5, UnitRipTable_L.(rsuf)(3) = UnitRipTable_L.(rsuf)(3)+1;
            case 2, UnitRipTable_L.(rsuf)(4) = UnitRipTable_L.(rsuf)(4)+1;
            case 3, UnitRipTable_L.(rsuf)(5) = UnitRipTable_L.(rsuf)(5)+1;
            case 4, UnitRipTable_L.(rsuf)(6) = UnitRipTable_L.(rsuf)(6)+1;
        end

          switch thisUnits.Selectivity_RScene(uid)
            case 0, UnitRipTable_R.(rsuf)(1) = UnitRipTable_R.(rsuf)(1)+1;
            case 1, UnitRipTable_R.(rsuf)(2) = UnitRipTable_R.(rsuf)(2)+1;
            case 0.5, UnitRipTable_R.(rsuf)(3) = UnitRipTable_R.(rsuf)(3)+1;
            case 2, UnitRipTable_R.(rsuf)(4) = UnitRipTable_R.(rsuf)(4)+1;
            case 3, UnitRipTable_R.(rsuf)(5) = UnitRipTable_R.(rsuf)(5)+1;
            case 4, UnitRipTable_R.(rsuf)(6) = UnitRipTable_R.(rsuf)(6)+1;
        end


          switch thisUnits.Selectivity_LScene(uid)
            case 0, UnitRipTable_C.(rsuf)(1) = UnitRipTable_C.(rsuf)(1)+1;
            case 1, UnitRipTable_C.(rsuf)(2) = UnitRipTable_C.(rsuf)(2)+1;
            case 0.5, UnitRipTable_C.(rsuf)(3) = UnitRipTable_C.(rsuf)(3)+1;
            case 2, UnitRipTable_C.(rsuf)(4) = UnitRipTable_C.(rsuf)(4)+1;
            case 3, UnitRipTable_C.(rsuf)(5) = UnitRipTable_C.(rsuf)(5)+1;
            case 4, UnitRipTable_C.(rsuf)(6) = UnitRipTable_C.(rsuf)(6)+1;
          end
    end

end

UnitRipTable_L.Sp_NSp = UnitRipTable_L.Sp_NSp_L+UnitRipTable_L.Sp_NSp_R+UnitRipTable_L.Sp_NSp_C+...
    UnitRipTable_L.Sp_NSp_L_R+UnitRipTable_L.Sp_NSp_L_C+UnitRipTable_L.Sp_NSp_R_C+UnitRipTable_L.Sp_NSp_L_R_C;
UnitRipTable_L.X_NSp = UnitRipTable_L.X_NSp_L+UnitRipTable_L.X_NSp_R+UnitRipTable_L.X_NSp_C+...
    UnitRipTable_L.X_NSp_L_R+UnitRipTable_L.X_NSp_L_C+UnitRipTable_L.X_NSp_R_C+UnitRipTable_L.X_NSp_L_R_C;

UnitRipTable_R.Sp_NSp = UnitRipTable_R.Sp_NSp_L+UnitRipTable_R.Sp_NSp_R+UnitRipTable_R.Sp_NSp_C+...
    UnitRipTable_R.Sp_NSp_L_R+UnitRipTable_R.Sp_NSp_L_C+UnitRipTable_R.Sp_NSp_R_C+UnitRipTable_R.Sp_NSp_L_R_C;
UnitRipTable_R.X_NSp = UnitRipTable_R.X_NSp_L+UnitRipTable_R.X_NSp_R+UnitRipTable_R.X_NSp_C+...
    UnitRipTable_R.X_NSp_L_R+UnitRipTable_R.X_NSp_L_C+UnitRipTable_R.X_NSp_R_C+UnitRipTable_R.X_NSp_L_R_C;

UnitRipTable_C.Sp_NSp = UnitRipTable_C.Sp_NSp_L+UnitRipTable_C.Sp_NSp_R+UnitRipTable_C.Sp_NSp_C+...
    UnitRipTable_C.Sp_NSp_L_R+UnitRipTable_C.Sp_NSp_L_C+UnitRipTable_C.Sp_NSp_R_C+UnitRipTable_C.Sp_NSp_L_R_C;
UnitRipTable_C.X_NSp = UnitRipTable_C.X_NSp_L+UnitRipTable_C.X_NSp_R+UnitRipTable_C.X_NSp_C+...
    UnitRipTable_C.X_NSp_L_R+UnitRipTable_C.X_NSp_L_C+UnitRipTable_C.X_NSp_R_C+UnitRipTable_C.X_NSp_L_R_C;
