Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Mother '\Processed Data\ripples_mat\R4'];
ROOT.Fig = [ROOT.Mother '\Processed Data\units_mat\ProfilingSheet\U3'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];

dir = '';

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


%%
thisRegion = 'SUB';
thisRegion2 = 'SUB_field';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_forAnalysis_RDI.xlsx']);
% UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion2 '.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM'};
unit = 5.0000e-04; ti = 1;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.05*Params.Fs;
dur = 0.4*Params.Fs;
thisFRMapSCALE=2;
Params.tbinDuration = 0.005;
filter_ns = 'C';


%%

%%
for uid=1:size(UnitsTable_A,1)
    try
        thisUnit = UnitsTable_A(uid,:);

        UnitID = thisUnit.ID{1};
        [thisRID,thisSID,thisTTID,thisCLID, thisFLID] = parsing_clusterID(UnitID,2);
        [thisRIDn,thisSIDn,thisTTIDn,thisCLIDn, thisFLIDn] = parsing_clusterID(UnitID,8);
        thisRSID = [thisRID '-' thisSID];

        if ~strcmp(thisRSID, thisRSID_old)
            Pos = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\ParsedPosition.mat']);
            diverging_point = get_divergingPoint(ROOT.Info, thisRID, thisSID);

            diverging_point = diverging_point*0.23;
            stem_end_index =  (max(Pos.y)-diverging_point)/thisFRMapSCALE;

            Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);
            Spike_f=LoadSpikeData_field(ROOT, thisRSID, [1:24], Params.cellfindn);
            disp([thisRSID ' plotting...'])
            thisRSID_old = thisRSID;
        end
        %% Load units

        Clist=jet(256);
        Clist_f = {'#FF6A68','#00B455','#00A2FF','#FA7344','#00B7A1'};

        Spk = Spike.(['TT' num2str(thisTTIDn)]).(['Unit' num2str(thisCLIDn)]);
        Maps = load([ROOT.Raw.Map '\rat' [thisRSID '-' num2str(thisTTIDn) '-' thisCLID] '.mat']);


        %% Load fields
        fids = find(UnitsTable_B.rat==thisRIDn & UnitsTable_B.session==thisSIDn & UnitsTable_B.TT==thisTTIDn & UnitsTable_B.AvgPeaktoValley==thisUnit.AvgPeaktoValley);

        clear Spk_f Maps_f
        for f = 1:size(fids,1)
            thisField = UnitsTable_B(fids(f),:);
            [thisRIDf,thisSIDf,thisTTIDf,thisCLIDf, thisFLIDf] = parsing_clusterID(thisField.ID{1},2);

            Spk_f(f) = Spike_f.(['TT' num2str(str2double((thisTTIDf)))]).(['Unit' num2str(str2double(thisCLIDf)) '_' num2str(str2double(thisFLIDf))]);
            Maps_f(f) = load([ROOT.Raw.Map '\rat' [thisRIDf '-' thisSIDf '-' num2str(str2double((thisTTIDf))) '-' thisCLIDf '-' thisFLIDf] '.mat']);

        end

        mx=max(Maps.skaggsMap1D{1});


        %%
        figure('position',[317,63,1400,915],'color','w');

        %% title
        subplot(3,2,1)
        title([UnitID ', ' thisUnit.experimenter{1}],'fontsize',15)

                text(-1,mx,[{'peak FR'} {'mean FR'} {'SI'}])
          text(5,mx,[{jjnum2str(Maps.onmazeMaxFR1D(1),2)} '' {jjnum2str(Maps.onmazeAvgFR1D(1),2)}...
              {jjnum2str(Maps.SpaInfoScore1D(1),2)}],'color','k')
       
        for f = 1:size(fids,1)
            text((f+1)*5,mx,[{jjnum2str(Maps_f(f).onmazeMaxFR1D(1),2)} '' {jjnum2str(Maps_f(f).onmazeAvgFR1D(1),2)}...
                {jjnum2str(Maps_f(f).SpaInfoScore1D(1),2)}],'color',hex2rgb(Clist_f{f}))
        end
xlim([0 30]); ylim([0 mx*2])

        axis off
        %% overall
        subplot(3,2,3)
        hold on
        plot(Maps.skaggsMap1D{1},'linewidth',4,'color','k')
        for f = 1:size(fids,1)
            plot(Maps_f(f).skaggsMap1D{1},'linewidth',4,'color',hex2rgb(Clist_f{f}))
        end

        settick(stem_end_index,max(Maps.skaggsMap1D{1}))
        title('overall')
        %% left
        subplot(3,2,2)
        hold on
        text(1,mx,[{'RDI'} 'p'],'color','k')
        %            plot(Maps.skaggsMap1D{2},'linewidth',4,'color','k')
        %            plot(Maps.skaggsMap1D{3},'linewidth',4,'color',[.7 .7 .7])
        %            plot(Maps.skaggsMap1D{3},'linewidth',4,'color','k','linestyle',':')
        for f = 1:size(fids,1)
            plot(Maps_f(f).skaggsMap1D{2},'linewidth',4,'color',hex2rgb(Clist_f{f}))
            p = plot(Maps_f(f).skaggsMap1D{3},'linewidth',4,'color',hex2rgb(Clist_f{f})); p.Color(4)=0.2;
            plot(Maps_f(f).skaggsMap1D{3},'linewidth',4,'color',hex2rgb(Clist_f{f}),'linestyle',':')
            text(f*5,mx,[{jjnum2str(UnitsTable_B.RDI_LScene(fids(f)),2)} '' jjnum2str(UnitsTable_B.RateP_LScene(fids(f)),2)],'color',hex2rgb(Clist_f{f}))
        end

        settick(stem_end_index,mx)
        title('Left scene (Zebra vs. Bamboo)')
        %% right
        subplot(3,2,4)
        hold on
        text(1,mx,[{'RDI'} 'p'],'color','k')

        for f = 1:size(fids,1)
            plot(Maps_f(f).skaggsMap1D{4},'linewidth',4,'color',hex2rgb(Clist_f{f}))
            p = plot(Maps_f(f).skaggsMap1D{5},'linewidth',4,'color',hex2rgb(Clist_f{f})); p.Color(4)=0.2;
            plot(Maps_f(f).skaggsMap1D{5},'linewidth',4,'color',hex2rgb(Clist_f{f}),'linestyle',':')
            text(f*5,mx,[{jjnum2str(UnitsTable_B.RDI_RScene(fids(f)),2)} '' jjnum2str(UnitsTable_B.RateP_RScene(fids(f)),2)],'color',hex2rgb(Clist_f{f}))
        end

        settick(stem_end_index,mx)
        title('Right scene (Pebbles vs. Mountains)')
        %% choice
        subplot(3,2,6)
        rectangle('position',[stem_end_index 0 48-stem_end_index mx], 'facecolor',[.8 .8 .8],'edgecolor','w')

        hold on
        text(1,mx,[{'RDI'} 'p'],'color','k')

        for f = 1:size(fids,1)
            plot(Maps_f(f).skaggsMap_left1D,'linewidth',4,'color',hex2rgb(Clist_f{f}))
            p = plot(Maps_f(f).skaggsMap_right1D,'linewidth',4,'color',hex2rgb(Clist_f{f})); p.Color(4)=0.2;
            plot(Maps_f(f).skaggsMap_right1D,'linewidth',4,'color',hex2rgb(Clist_f{f}),'linestyle',':')
            text(f*5,mx,[{jjnum2str(UnitsTable_B.RDI_LR(fids(f)),2)} '' jjnum2str(UnitsTable_B.RateP_LR(fids(f)),2)],'color',hex2rgb(Clist_f{f}))
        end

        settick(stem_end_index,mx)
        title('Choice (Left vs. Right)')

        %% save fig
       if ~exist([ROOT.Fig '\'  thisRegion]), mkdir([ROOT.Fig '\'  thisRegion]); end
        saveas(gca,[ROOT.Fig '\'  thisRegion '\' UnitID '.svg'])
        saveas(gca,[ROOT.Fig '\'  thisRegion '\' UnitID '.png'])
        close all

    catch
        close all
    end
end


function settick(Dv,mx)
xlabel('position')
xticks([0 Dv 48])
xticklabels({'StBox','Dv','Fd'})
xlim([0 48])
ylabel('firing rate (Hz)')
ylim([0 mx*2])
line([Dv Dv],[0 mx*2],'linestyle','--','color','k')
end