Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R5'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R6'];
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
RipplesTable_p = RipplesTable(RipplesTable.nPCs>2,:);
for r=1:size(RipplesTable_p,1)
    try
        thisRip = RipplesTable_p(r,:);
        thisReact = ReactTable(strcmp(ReactTable.RippleID,thisRip.ID{1}),:);
        thisUnits=table;

        [C,ia] = unique(thisReact.UnitID);
thisUnits = thisReact(ia,:);
thisUnits = sortrows(thisUnits,{'SpkTime_fromRippleStart'});

        temp_r = decimalToBinaryVector(thisRip.TTs);
        RippleTT = flip(length((temp_r))+1 - find(temp_r))';
        Exper = cell2mat(thisRip.experimenter);
        if ~ismember(Exper,Experimenter), continue; end
        thisRSID = [jmnum2str(thisRip.rat,3) '-' jmnum2str(thisRip.session,2)];
        load([ROOT.Rip ['\R-' thisRip.ID{1} '.mat']])

        clusters=imread([ROOT.Fig4 '\' thisRip.ID{1} '.png']);
        replay=imread([ROOT.Fig5 '\' thisRip.ID{1} '.png']);
        %
        figure('position',[100 100 2100 900],'color','w')
        subplot('position',[0 0 .45 1])
        imshow(clusters)
        %
        subplot('position',[.79 .02 .25 .8])
        imshow(replay(120:end,120:445,:))
        %

        subplot('position', [.45 .5 .15 .14])
        scatter(thisUnits.RDI_LScene,[1:size(thisUnits,1)],40,'k','filled')
        xlim([-1.5 1.5]); ylim([.5 size(thisUnits,1)+.5])
        line([0 0],[0 size(thisUnits,1)+1],'color','k')
        line([RDI_L.act_mean RDI_L.act_mean],[0 size(thisUnits,1)+1],'color','r')
        line([RDI_L.act_median RDI_L.act_median],[0 size(thisUnits,1)+1],'color','b')
        title('Left scene selectivity')
        axis ij

        subplot('position', [.45 .3 .15 .14])
                scatter(thisUnits.RDI_RScene,[1:size(thisUnits,1)],40,'k','filled')
        xlim([-1.5 1.5]); ylim([.5 size(thisUnits,1)+.5])
        line([0 0],[0 size(thisUnits,1)+1],'color','k')
        line([RDI_R.act_mean RDI_R.act_mean],[0 size(thisUnits,1)+1],'color','r')
        line([RDI_R.act_median RDI_R.act_median],[0 size(thisUnits,1)+1],'color','b')
        title('Right scene selectivity')
        axis ij

        subplot('position', [.45 .1 .15 .14])
                scatter(thisUnits.RDI_LR,[1:size(thisUnits,1)],40,'k','filled')
        xlim([-1.5 1.5]); ylim([.5 size(thisUnits,1)+.5])
        line([0 0],[0 size(thisUnits,1)+1],'color','k')
        line([RDI_C.act_mean RDI_C.act_mean],[0 size(thisUnits,1)+1],'color','r')
        line([RDI_C.act_median RDI_C.act_median],[0 size(thisUnits,1)+1],'color','b')
        title('Choice selectivity')
        axis ij
        %

        subplot('position', [.63 .66 .06 .01])
        title('mean RDI','FontSize',15); axis off
        subplot('position', [.72 .66 .06 .01])
        title('median RDI','FontSize',15); axis off
        
        subplot('position', [.63 .5 .06 .15])
        perm_hist(RDI_L,'mean')


        subplot('position', [.63 .3 .06 .15])
        perm_hist(RDI_R,'mean')

        subplot('position', [.63 .1 .06 .15])
        perm_hist(RDI_C,'mean')
        
        subplot('position', [.72 .5 .06 .15])
        perm_hist(RDI_L,'median')
        subplot('position', [.72 .3 .06 .15])
        perm_hist(RDI_R,'median')
        subplot('position', [.72 .1 .06 .15])
        perm_hist(RDI_C,'median')

        saveas(gca,[ROOT.Fig '\' thisRip.ID{1} '.png'])
        close all


    catch
        close all
    end

end

function perm_hist(RDI,perm)
if RDI.(['p_' perm])<0.05 && ~isnan(RDI.(['act_' perm])), cl='r'; else, cl='k'; end
if strcmp(perm,'mean'), m='m'; cl2='r'; else, m='M'; cl2='b'; end

if ~isnan(RDI.(['act_' perm]))
histogram(RDI.dist.(perm),'facecolor','k','BinWidth', .05)
line([RDI.(['act_' perm]) RDI.(['act_' perm])], [0 1000],'color',cl2)
end

title(['p= ' jjnum2str(RDI.(['p_' perm]),3) ', ' m '= ' jjnum2str(RDI.(['act_' perm]),3)],'color',cl)
xlim([-.7 .7])
end