Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R5'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R8'];
ROOT.Fig4 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R4'];
ROOT.Fig5 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R5'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end
if ~exist(ROOT.Rip), mkdir(ROOT.Rip); end

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
randN=2000;
%%
RipplesTable_p = RipplesTable(RipplesTable.nPCs>4,:);


for r=1:size(RipplesTable_p,1)
    try
        thisRip = RipplesTable_p(r,:);
        thisReact = ReactTable(strcmp(ReactTable.RippleID,thisRip.ID{1}),:);
        thisUnits=table;



        %         [C,ia] = unique(thisReact.UnitID);
        % thisUnit = thisReact(ia,:);
        % thisUnit = sortrows(thisUnit,{'UnitID'});
        % [~,ord] = sortrows(thisUnit,{'SpkTime_fromRippleStart'});

        temp_r = decimalToBinaryVector(thisRip.TTs);
        RippleTT = flip(length((temp_r))+1 - find(temp_r))';
        Exper = cell2mat(thisRip.experimenter);
        if ~ismember(Exper,Experimenter), continue; end
        thisRSID = [jmnum2str(thisRip.rat,3) '-' jmnum2str(thisRip.session,2)];
        load([ROOT.Rip ['\R-' thisRip.ID{1} '.mat']])
        thisPool = UnitsTable(UnitsTable.rat==RipplesTable_p.rat(r) & UnitsTable.session==RipplesTable_p.session(r),:);


        RipplesTable_p.nRDIs_L(r)=sum(~isnan(thisUnits.RDI_LScene));
        RipplesTable_p.nRDIs_R(r)=sum(~isnan(thisUnits.RDI_RScene));
        RipplesTable_p.nRDIs_C(r)=sum(~isnan(thisUnits.RDI_LR));
        RipplesTable_p.nRDIsMax(r) = nanmax([RipplesTable_p.nRDIs_L(r),RipplesTable_p.nRDIs_R(r),RipplesTable_p.nRDIs_C(r)]);
        RipplesTable_p.nRDIsMin(r) = nanmin([RipplesTable_p.nRDIs_L(r),RipplesTable_p.nRDIs_R(r),RipplesTable_p.nRDIs_C(r)]);
        if RipplesTable_p.nRDIsMax(r)<5, continue;end

        clusters=imread([ROOT.Fig4 '\' thisRip.ID{1} '.png']);
        replay=imread([ROOT.Fig5 '\' thisRip.ID{1} '.png']);
        %
        figure('position',[100 100 2100 900],'color','w')

        subplot('position',[.05 .95 .14 .05])
        text(0,0,thisRip.ID,'fontsize',15,'fontweight','b')
        axis off
        subplot('position',[.9 .95 .1 .05])
        text(0, 0, ['printed on ' date],'FontSize',12);
        axis off

        x1=.63;
        w1=.81;
        w2=.9;

        subplot('position',[.15 .02 .45 .9])
        imshow(clusters(70:end-100,100:end,:))
        %
        subplot('position',[0 .02 .14 .7])
        imshow(replay(120:end,120:430,:))
        %

        subplot('position', [x1 .5 .15 .14]); hold on
        RDI_L = scatter_RDI(RDI_L,thisUnits,thisPool,'RDI_LScene',20);

        subplot('position', [x1 .3 .15 .14]); hold on
        RDI_R = scatter_RDI(RDI_R, thisUnits,thisPool,'RDI_RScene',20);


        subplot('position', [x1 .1 .15 .14]); hold on
        RDI_C = scatter_RDI(RDI_C,thisUnits,thisPool,'RDI_LR',20);
        %

        subplot('position', [w1 .66 .06 .01])
        title('mean RDI','FontSize',15); axis off
        subplot('position', [w2 .66 .06 .01])
        title('median RDI','FontSize',15); axis off

        subplot('position', [w1 .5 .06 .15])
        RDI_L = perm_hist(RDI_L,'mean');


        subplot('position', [w1 .3 .06 .15])
        RDI_R = perm_hist(RDI_R,'mean');

        subplot('position', [w1 .1 .06 .15])
        RDI_C = perm_hist(RDI_C,'mean');

        subplot('position', [w2 .5 .06 .15])
        RDI_L = perm_hist(RDI_L,'median');
        subplot('position', [w2 .3 .06 .15])
        RDI_R = perm_hist(RDI_R,'median');
        subplot('position', [w2 .1 .06 .15])
        RDI_C = perm_hist(RDI_C,'median');
                RipplesTable_p.M_m(r) = max([abs(RDI_L.act_mean-RDI_L.act_median),abs(RDI_R.act_mean-RDI_R.act_median),abs(RDI_C.act_mean-RDI_C.act_median)]);



        ROOT.Fig_en = [ROOT.Fig '\Ensemble_' num2str(thisRip.nRDIsMax)];
        if ~exist(ROOT.Fig_en), mkdir(ROOT.Fig_en); end
        saveas(gca,[ROOT.Fig_en '\' thisRip.ID{1} '.png'])
        close all


        RipplesTable_p.pRDI_L(r) = RDI_L.p_mean;
        RipplesTable_p.pRDI_R(r) = RDI_R.p_mean;
        RipplesTable_p.pRDI_C(r) = RDI_C.p_mean;
    catch
        close all
    end

end

RipplesTable_c = RipplesTable_p(RipplesTable_p.nRDIsMax>=5,:);
writetable(RipplesTable_c,[ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx'])

RipplesTable_c = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);

function RDI = perm_hist(RDI,perm)
RDI.dist.(perm)=RDI.dist.(perm)(~isnan(RDI.dist.(perm)));
RDI.(['p_' perm]) = min(sum(RDI.(['act_' perm])<RDI.dist.(perm)),sum(RDI.(['act_' perm])>RDI.dist.(perm)))/length(RDI.dist.(perm));

if RDI.(['p_' perm])<0.05 && ~isnan(RDI.(['act_' perm])), cl='r'; else, cl='k'; end
if strcmp(perm,'mean'), m='m'; cl2='r'; else, m='M'; cl2='b'; end

if ~isnan(RDI.(['act_' perm]))
    histogram(RDI.dist.(perm),'facecolor','k','BinWidth', .05)
    line([RDI.(['act_' perm]) RDI.(['act_' perm])], [0 1000],'color',cl2)
end

title(['p= ' jjnum2str(RDI.(['p_' perm]),3) ', ' m '= ' jjnum2str(RDI.(['act_' perm]),3)],'color',cl)
xlim([-.7 .7])
end

function RDI = scatter_RDI(RDI,thisUnits,thisPool,field,n)
switch field
    case 'RDI_LScene', titp = 'Left scene';
    case 'RDI_RScene', titp = 'Right scene';
    case 'RDI_LR', titp = 'Choice';
end

PoolL = sort(thisPool.(field)); PoolL=PoolL(~isnan(PoolL));
thisUnits = sortrows(thisUnits,{field});
id = (thisUnits.(field)<nanmean(thisUnits.(field))+n*nanstd(thisUnits.(field))) &...
    (thisUnits.(field)>nanmean(thisUnits.(field))-n*nanstd(thisUnits.(field)));
thisUnits_in = thisUnits(id,:);
thisUnits_out = thisUnits(~id,:);
scatter(PoolL,[1:size(PoolL,1)],40,'k','filled','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
scatter(thisUnits_in.(field),knnsearch(PoolL,thisUnits_in.(field)),40,'r','filled')
scatter(thisUnits_out.(field),knnsearch(PoolL,thisUnits_out.(field)),40,'r')
xlim([-1.5 1.5]); ylim([.5,size(PoolL,1)+.5])
line([0 0],[.5,size(PoolL,1)+.5],'color','k', 'linestyle','--')
line([nanmean(thisUnits_in.(field)) nanmean(thisUnits_in.(field))],[.5,size(PoolL,1)+.5],'color','r')
line([nanmedian(thisUnits_in.(field)) nanmedian(thisUnits_in.(field))],[.5,size(PoolL,1)+.5],'color','b')
title([titp ' selectivity stdev: ' jjnum2str(nanstd(thisUnits.(field)),3) '/'...
    jjnum2str(nanstd(thisPool.(field)),3) '=' jjnum2str(nanstd(thisUnits.(field))/nanstd(thisPool.(field)),3)])
title([titp ' selectivity distribution'])
axis ij

RDI.act_mean = nanmean(thisUnits_in.(field));
RDI.act_median = nanmedian(thisUnits_in.(field));



end
