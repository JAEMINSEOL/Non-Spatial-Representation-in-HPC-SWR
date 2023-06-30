thisRegion = 'SUB';
thisRegion2 = 'SUB_field';
RipplesTable_sub = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_RDIs_UV_cell_HeteroIn_AllPopul' '.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);

UnitsTable_sub_f = readtable([ROOT.Units '\UnitsTable_' 'SUB_field' '_forAnalysis.xlsx']);
UnitsTable_ca1_f = readtable([ROOT.Units '\UnitsTable_' 'CA1_field' '_forAnalysis.xlsx']);


UnitsTable_sub = readtable([ROOT.Units '\UnitsTable_' 'SUB' '_forAnalysis.xlsx']);
UnitsTable_ca1 = readtable([ROOT.Units '\UnitsTable_' 'CA1' '_forAnalysis.xlsx']);
% writetable(RipplesTable,[ROOT.Save '\RipplesTable_' thisRegion2 '_RDIs_UV_cell' '.xlsx'],'writemode','replacefile')

CList = [[139 193 69]/255; [29 111 169]/255];
%%
for rid=1:size(RipplesTable,1)
    thisRip =RipplesTable(rid,:);
    if thisRip.DecodingP_all<0.05
        rsuf1 = 'Sp';
    else
        rsuf1 = 'X';
    end
    if nanmin([thisRip.pRDI_L_UV,thisRip.pRDI_R_UV,thisRip.pRDI_C_UV])<0.05
        rsuf2 = 'NSp';
        if thisRip.pRDI_L_UV<0.05, rsufl='_L'; else, rsufl = ''; end
        if thisRip.pRDI_R_UV<0.05, rsufr='_R'; else, rsufr = ''; end
        if thisRip.pRDI_C_UV<0.05, rsufc='_C'; else, rsufc = ''; end
    else
        rsuf2 = 'X'; rsufl = ''; rsufr = ''; rsufc = '';
    end

    rsuf = [rsuf1 '_' rsuf2 rsufl rsufr rsufc];

    RipplesTable.Type{rid} = rsuf;
end

%%

thisRSID_p='';
sid=0;
SList=[];

for rid=1:size(RipplesTable,1)
     thisRip =RipplesTable(rid,:);
    rsuf = thisRip.Type{1};
    thisRSID = [jmnum2str(thisRip.rat,3) '-' jmnum2str(thisRip.session,2)];
    if ~strcmp(thisRSID, thisRSID_p)
        sid=sid+1;
        SList=[SList;thisRSID];
        thisRSID_p = thisRSID;
    end
end

SessionRipType = table('Size',[size(SList,1) 17], 'VariableTypes', {'string' 'double' 'double' 'double' 'double',...
    'double' 'double' 'double' 'double' 'double' 'double' 'double' 'double' 'double' 'double' 'double' 'double'},...
    'VariableNames',{'Session','X_X','Sp_X',...
    'Sp_NSp_L' 'Sp_NSp_R' 'Sp_NSp_C' 'Sp_NSp_L_R' 'Sp_NSp_L_C' 'Sp_NSp_R_C' 'Sp_NSp_L_R_C',...
    'X_NSp_L' 'X_NSp_R' 'X_NSp_C' 'X_NSp_L_R' 'X_NSp_L_C' 'X_NSp_R_C' 'X_NSp_L_R_C'});
SessionRipType.Session = SList;
thisRSID_p='';
sid=0;
for rid=1:size(RipplesTable,1)
     thisRip =RipplesTable(rid,:);
    rsuf = thisRip.Type{1};
   thisRSID = [jmnum2str(thisRip.rat,3) '-' jmnum2str(thisRip.session,2)];
    if ~strcmp(thisRSID, thisRSID_p)
        sid=sid+1;
        thisRSID_p = thisRSID;
    end

    SessionRipType.(rsuf)(sid) = SessionRipType.(rsuf)(sid)+1;

end

SessionRipType.Sp_NSp = SessionRipType.Sp_NSp_L+SessionRipType.Sp_NSp_R+SessionRipType.Sp_NSp_C+...
    SessionRipType.Sp_NSp_L_R+SessionRipType.Sp_NSp_L_C+SessionRipType.Sp_NSp_R_C+SessionRipType.Sp_NSp_L_R_C;
SessionRipType.X_NSp = SessionRipType.X_NSp_L+SessionRipType.X_NSp_R+SessionRipType.X_NSp_C+...
    SessionRipType.X_NSp_L_R+SessionRipType.X_NSp_L_C+SessionRipType.X_NSp_R_C+SessionRipType.X_NSp_L_R_C;

%%
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

%%
RipPartL =table; RipPartR = table; RipPartC = table;
for rid=1:size(RipplesTable,1)
     thisRip =RipplesTable(rid,:);
         thisReact_A = ReactTable(strcmp(ReactTable.RippleID,thisRip.ID{1}),:);
           thisUnits=table;
    for u=1:size(thisReact_A,1)
        thisUnits =[thisUnits;UnitsTable_A(find(cellfun(Params.cellfind(thisReact_A.UnitID(u)),UnitsTable_A.ID)),:)];
    end
    rsuf = thisRip.Type{1};

    RipPartL.SF_N(rid) = sum(thisUnits.Selectivity_LScene==0)/size(thisUnits,1);
    RipPartL.SF_S(rid) = sum(thisUnits.Selectivity_LScene==1)/size(thisUnits,1);
    RipPartL.MF_N(rid) = sum(thisUnits.Selectivity_LScene==0.5)/size(thisUnits,1);
    RipPartL.MF_S(rid) = sum(thisUnits.Selectivity_LScene==2)/size(thisUnits,1);
    RipPartL.MF_Homo(rid) = sum(thisUnits.Selectivity_LScene==3)/size(thisUnits,1);
    RipPartL.MF_Hetero(rid) = sum(thisUnits.Selectivity_LScene==4)/size(thisUnits,1);

    RipPartR.SF_N(rid) = sum(thisUnits.Selectivity_RScene==0)/size(thisUnits,1);
    RipPartR.SF_S(rid) = sum(thisUnits.Selectivity_RScene==1)/size(thisUnits,1);
    RipPartR.MF_N(rid) = sum(thisUnits.Selectivity_RScene==0.5)/size(thisUnits,1);
    RipPartR.MF_S(rid) = sum(thisUnits.Selectivity_RScene==2)/size(thisUnits,1);
    RipPartR.MF_Homo(rid) = sum(thisUnits.Selectivity_RScene==3)/size(thisUnits,1);
    RipPartR.MF_Hetero(rid) = sum(thisUnits.Selectivity_RScene==4)/size(thisUnits,1);

    RipPartC.SF_N(rid) = sum(thisUnits.Selectivity_LR==0)/size(thisUnits,1);
    RipPartC.SF_S(rid) = sum(thisUnits.Selectivity_LR==1)/size(thisUnits,1);
    RipPartC.MF_N(rid) = sum(thisUnits.Selectivity_LR==0.5)/size(thisUnits,1);
    RipPartC.MF_S(rid) = sum(thisUnits.Selectivity_LR==2)/size(thisUnits,1);
    RipPartC.MF_Homo(rid) = sum(thisUnits.Selectivity_LR==3)/size(thisUnits,1);
    RipPartC.MF_Hetero(rid) = sum(thisUnits.Selectivity_LR==4)/size(thisUnits,1);
         
end

 RipPartL.RipType = RipplesTable.Type;
 RipPartR.RipType = RipplesTable.Type;
 RipPartC.RipType = RipplesTable.Type;
%%
thisRegion = 'CA1';
Ca1 = load([ROOT.Save '\RipPart_' thisRegion '.mat'], 'RipPartL','RipPartR','RipPartC');
thisRegion = 'SUB';
Sub=load([ROOT.Save '\RipPart_' thisRegion '.mat'], 'RipPartL','RipPartR','RipPartC');

T1=Sub.RipPartC;
T2=Ca1.RipPartC;
% 
% T1 = T1(find(cellfun(Params.cellfindn2('X_NSp'),table2array(T1(:,7)))),:);
% T2 = T2(find(cellfun(Params.cellfindn2('X_NSp'),table2array(T2(:,7)))),:);


figure;
data = [mean(T1.SF_N) mean(T1.MF_N) mean(T1.SF_S) mean(T1.MF_S) mean(T1.MF_Homo) mean(T1.MF_Hetero);...
    mean(T2.SF_N) mean(T2.MF_N) mean(T2.SF_S) mean(T2.MF_S) mean(T2.MF_Homo) mean(T2.MF_Hetero)];

err = [std(T1.SF_N)/sqrt(size(T1.SF_N,1)) std(T1.MF_N)/sqrt(size(T1.MF_N,1)) std(T1.SF_S)/sqrt(size(T1.SF_S,1)),...
    std(T1.MF_S)/sqrt(size(T1.MF_S,1)) std(T1.MF_Homo)/sqrt(size(T1.MF_Homo,1)) std(T1.MF_Hetero)/sqrt(size(T1.MF_Hetero,1));...
    std(T2.SF_N)/sqrt(size(T2.SF_N,1)) std(T2.MF_N)/sqrt(size(T2.MF_N,1)) std(T2.SF_S)/sqrt(size(T2.SF_S,1)),...
    std(T2.MF_S)/sqrt(size(T2.MF_S,1)) std(T2.MF_Homo)/sqrt(size(T2.MF_Homo,1)) std(T2.MF_Hetero)/sqrt(size(T2.MF_Hetero,1))];

b = bar([1 2 3 4 5 6],data);
b(1).FaceColor = CList(1,:);
b(2).FaceColor = CList(2,:);
hold on
for k = 1:size(data,1)
    % get x positions per group
    xpos = b(k).XData + b(k).XOffset;
    % draw errorbar
    errorbar(xpos, data(k,:), err(k,:), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
end

t1 = UnitsTable_sub.Selectivity_LR;
t2 = UnitsTable_ca1.Selectivity_LR;
type_CA1 = [sum(t2==0) sum(t2==0.5) sum(t2==1) sum(t2==2) sum(t2==3) sum(t2==4)] / size(t2,1);
type_SUB = [sum(t1==0) sum(t1==0.5) sum(t1==1) sum(t1==2) sum(t1==3) sum(t1==4)] / size(t1,1);

% c = bar([1 2 3 4 5 6],[type_SUB;type_CA1]);
% c(1).FaceAlpha=0; c(1).EdgeColor='r';
% c(2).FaceAlpha=0; c(2).EdgeColor='r';

set(gca,'fontsize',12,'fontweight','b')

xticklabels({'SF-N' 'MF-N' 'SF-Selective' 'MF-Single' 'MF-Homo' 'MF-Hetero'})
legend({'SUB','CA1'});
ylabel('Proportion')
ylim([0 .6])

[~,p]  = ttest(T1.MF_N-type_SUB(2))
%%
sub =sum(table2array(SRT_SUB(:,2:19)),2); sub(sub<25)=nan;
ca1 =sum(table2array(SRT_CA1(:,2:19)),2); ca1(ca1<25)=nan;

[p,h]  = ranksum(SRT_CA1.Sp_X./ca1,SRT_SUB.Sp_X./sub)
[p,h]  = ranksum(SRT_CA1.X_NSp./ca1,SRT_SUB.X_NSp./sub)

%%

x1=[];x2=[];y1=[];y2=[];
cont = 'LScene';
T1 = UnitsTable_ca1;
x1 = T1.(['RDI_' cont])(T1.(['Selectivity_' cont])>2 | T1.(['Selectivity_' cont])==0.5);
y1 = T1.(['RDI_' cont])(T1.(['Selectivity_' cont])==1 | T1.(['Selectivity_' cont])==0);

T1 = UnitsTable_sub;
x2 = T1.(['RDI_' cont])(T1.Selectivity_LScene>2 | T1.Selectivity_LScene==0.5);
y2 = T1.(['RDI_' cont])(T1.Selectivity_LScene==1 | T1.Selectivity_LScene==0);

% cont = 'RScene';
% T1 = UnitsTable_ca1;
% x1 = [x1;T1.(['RDI_' cont])(T1.(['Selectivity_' cont])>2 | T1.(['Selectivity_' cont])==0.5)];
% y1 = [y1;T1.(['RDI_' cont])(T1.(['Selectivity_' cont])==1 | T1.(['Selectivity_' cont])==0)];
% 
% T1 = UnitsTable_sub;
% x2 = [x2;T1.(['RDI_' cont])(T1.Selectivity_LScene>2 | T1.Selectivity_LScene==0.5)];
% y2 = [y2; T1.(['RDI_' cont])(T1.Selectivity_LScene==1 | T1.Selectivity_LScene==0)];
% 
% cont = 'LR';
% T1 = UnitsTable_ca1;
% x1 = [x1;T1.(['RDI_' cont])(T1.(['Selectivity_' cont])>2 | T1.(['Selectivity_' cont])==0.5)];
% y1 = [y1;T1.(['RDI_' cont])(T1.(['Selectivity_' cont])==1 | T1.(['Selectivity_' cont])==0)];
% 
% T1 = UnitsTable_sub;
% x2 = [x2;T1.(['RDI_' cont])(T1.Selectivity_LScene>2 | T1.Selectivity_LScene==0.5)];
% y2 = [y2; T1.(['RDI_' cont])(T1.Selectivity_LScene==1 | T1.Selectivity_LScene==0)];
figure;
histogram(x2,'Normalization','probability','binwidth',0.05)
hold on
histogram(y2,'Normalization','probability','binwidth',0.05)
xlabel('scene/choice selectivity')

[h,p] = kstest2(x2,y2)
[h,p] = kstest2(x1,y1)
%%
figure;
data = [nanmean(abs([x2])) nanmean(abs([y2]));nanmean(abs([x1])) nanmean(abs([y1]))];
err = [nanstd(abs(x2))/sqrt(length(x2)) nanstd(abs(y2))/sqrt(length(y2));...
    nanstd(abs(x1))/sqrt(length(x1)) nanstd(abs(y1))/sqrt(length(y1))];
hb = bar(data);
hold on
for k = 1:size(data,2)
    % get x positions per group
    xpos = hb(k).XData + hb(k).XOffset;
    % draw errorbar
    errorbar(xpos, data(:,k), err(:,k), 'LineStyle', 'none', ... 
        'Color', 'k', 'LineWidth', 1);
end

xticklabels({'SUB','CA1'})
set(gca,'fontsize',12,'fontweight','b')
legend({'MF','SF'});
% title(cont)

[h,p] = ttest2(abs([x1]),abs([y1]))
[h,p] = ttest2(abs([x2]),abs([y2]))


%%
% RipplesTable_sub = readtable([ROOT.Save '\RipplesTable_' 'SUB_field' '_RDIs_UV_cell' '.xlsx']);
% RipplesTable_ca1 = readtable([ROOT.Save '\RipplesTable_' 'CA1_field' '_RDIs_UV_cell' '.xlsx']);
figure;
histogram(RipplesTable_sub.nFields,'FaceColor',[139 193 69]/255,'Normalization','probability');
hold on
histogram(RipplesTable_ca1.nFields,'FaceColor',[29 111 169]/255,'Normalization','probability');

xlabel('# of fields for each reactivation')
ylabel('proportion')
legend({'SUB','CA1'});
% 
% hold on
% histogram(UnitsTable_ca1.RDI_LScene(UnitsTable_ca1.Selectivity_LScene==1),'Normalization','probability','binwidth',0.1)


%%

x=UnitsTable_sub.SI; y=UnitsTable_ca1.SI;

figure;
subplot(1,3,[1 2])
histogram(x,'Normalization','probability','binwidth',0.2,'FaceColor',CList(1,:))
hold on
histogram(y,'Normalization','probability','binwidth',0.2,'FaceColor',CList(2,:))
set(gca,'fontsize',12,'fontweight','b')
xlabel('SI score')
ylabel('proportion')
legend({'SUB','CA1'});

subplot(1,3,3)
b=bar([1 2], [mean(x) mean(y)]);
b.FaceColor='flat';
b.CData(1,:)=CList(1,:);
b.CData(2,:)=CList(2,:);
hold on
er = errorbar([1 2], [mean(x) mean(y)], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))], [std(x)/sqrt(length(x)) std(y)/sqrt(length(y))]);
er.Color = [0 0 0];                            
er.LineStyle = 'none';
set(gca,'fontsize',12,'fontweight','b')
xticklabels({'SUB','CA1'})
ylabel('SI score')

