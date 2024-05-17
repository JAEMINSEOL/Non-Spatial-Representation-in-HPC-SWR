Initial_SWRFilter_common;
warning off


RegionList = {'SUB','CA1'};
for reg=2:2
thisR = RegionList{reg};
% 

thisRegion0 = thisR;
thisRegion = thisR;
thisRegion2 = [thisR '_field'];

% thisRegion0 = 'CA1';
% thisRegion = 'CA1';
% thisRegion2 = 'CA1_field';

% thisRegion0 = 'SUB';
% thisRegion = 'SUB';
% thisRegion2 = 'SUB_field';


ROOT.Save = [ROOT.Processed ''];
ROOT.Rip0 = [ROOT.Processed '\ripples_mat\R0\' thisRegion];
ROOT.Rip = [ROOT.Processed '\ripples_mat\R2'];
ROOT.Rip4 = [ROOT.Processed '\ripples_mat\R4'];
ROOT.Rip5 = [ROOT.Processed '\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Processed '\ripples_mat\ProfilingSheet\R41 (raster only)_' thisRegion];
ROOT.Units = [ROOT.Processed '\units_mat\U2'];
ROOT.Behav = [ROOT.Processed '\behavior_mat'];

dir = '';

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


%%
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis_final.xlsx']);
RipplesTable2 = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis_240310.xlsx']);
% UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion2 '.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 1;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.05*Params.Fs;
dur = 0.4*Params.Fs;
thisFRMapSCALE=2;
Params.tbinDuration = 0.005;

filter_ns = 'C';
filter_ns2 = 'LR';
for sid=1:size(RipplesTable,1)
try
  thisRip = RipplesTable(sid,:);
        %                 if thisRip.correctness==1, continue; end
        thisRID = jmnum2str(thisRip.rat,3);
        thisSID = jmnum2str(thisRip.session,2);
        %         if ~strcmp(thisRID,'561'), continue; end



        Ist = thisRip.STindex-mar;
        Ied = thisRip.EDindex+mar;
        temp_r = decimalToBinaryVector(thisRip.TTs);
        RippleTT = flip(length((temp_r))+1 - find(temp_r))';
        Exper = cell2mat(thisRip.experimenter);
        if ~ismember(Exper,Experimenter), continue; end
        thisRSID = [jmnum2str(thisRip.rat,3) '-' jmnum2str(thisRip.session,2)];


        if ~strcmp(thisRSID, thisRSID_old)
            Pos = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\ParsedPosition.mat']);
            diverging_point = get_divergingPoint(ROOT.Info, thisRID, thisSID);

            diverging_point = diverging_point*0.23;
            stem_end_index =  (max(Pos.y)-diverging_point)/thisFRMapSCALE;

            Recording_region_TT = Recording_region(thisRSID,:);
            TargetTT = GetTargetTT(ROOT,thisRSID,thisRegion,Params,0.03);

            thisTT_table = TT_table(TT_table.rat==thisRip.rat & TT_table.session==thisRip.session,:);
            for t=1:size(thisTT_table,1)
                if ~ismember(thisTT_table.TT(t),TargetTT)
                    thisTT_table.TT(t)=0;
                end
            end
            thisTT_table= thisTT_table(thisTT_table.TT~=0,:);
            [~,t] = max(thisTT_table.RippleBandMean);
            TargetTT_p = thisTT_table.TT(t);

            EEG = LoadEEGData(ROOT, thisRSID, TargetTT_p,Params,Params_Ripple);


            clusters_A = UnitsTable_A(UnitsTable_A.rat==thisRip.rat & UnitsTable_A.session==thisRip.session,:);
            clusters_B = UnitsTable_B(UnitsTable_B.rat==thisRip.rat & UnitsTable_B.session==thisRip.session,:);
            Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);
            load([ROOT.Rip0 '\' thisRSID '.mat'])
            disp([thisRSID ' plotting...'])
            thisRSID_old = thisRSID;
            clusters_B.RDI_LR(clusters_B.PeakBin>stem_end_index)=nan;
            %% remove heterogeneous cells

            %              clusters_B = clusters_B(clusters_B.(['Selectivity_' filter_ns2])>=1 & clusters_B.(['Selectivity_' filter_ns2])<4,:);
            %              clusters_A = clusters_A(clusters_A.(['Selectivity_' filter_ns2])>=1 & clusters_A.(['Selectivity_' filter_ns2])<4,:);
            %% Load FRMap

        end

%%
if nanmin([thisRip.pBinomDev_L_UV, thisRip.pBinomDev_R_UV,thisRip.pBinomDev_C_UV])<0.05
    suf1=['TR'];
else
    suf1='X';
end

        [aspks_epoch,aspks_epoch_in,Unitsa,UnitsAa] = LoadUnits(thisRip,clusters_A,Spike,Params,mar);
clusters_A.isMF = clusters_A.NumField>1;
clusters_A.num=[1:size(clusters_A,1)]';
        c = sortrows(clusters_A,'isMF');
        c.num2 = [1:size(clusters_A,1)]';
        [clusters_A,ia] = sortrows(c,'num');

        aspks_epoch(:,11) = clusters_A.num2(aspks_epoch(:,2));
        %%
figure('position',[145,128,370,700],'color','w');

            %% EEG
        thisEEG = EEG.(['TT' num2str(TargetTT_p)]).Raw(Ist:Ied);
        thisEEG_Ripple = EEG.(['TT' num2str(TargetTT_p)]).Filtered(Ist:Ied);
        hold on

        subplot(5,1,1)
        plot(thisEEG,'k')
                        title(['(' thisRegion0 ') ' cell2mat(thisRip.ID) ', ' suf1],...
            'fontsize',15,'Interpreter','none')
        x1=mar; x2=mar+thisRip.RippleDuration*Params.Fs;
        dur2=x2+mar;
        xlim([0 dur2])
        line([x1 x1], [min(thisEEG) max(thisEEG)+50], 'color','r','linestyle','--')
        line([x2 x2], [min(thisEEG) max(thisEEG)+50], 'color','r','linestyle','--')


        subplot(5,1,2)
        plot(thisEEG_Ripple,'b')

        x1=mar; x2=mar+thisRip.RippleDuration*Params.Fs;
        dur2=x2+mar;
        xlim([0 dur2])
        line([x1 x1], [min(thisEEG_Ripple) max(thisEEG_Ripple)+50], 'color','r','linestyle','--')
        line([x2 x2], [min(thisEEG_Ripple) max(thisEEG_Ripple)+50], 'color','r','linestyle','--')


                %% reactivation rasterplot
        subplot(5,1,[3:5])
        hold on
        yl = size(clusters_A,1)+1;
        x1=mar/2; x2=thisRip.RippleDuration*1e3+mar/2;
        xlim([0 dur2/2])
        ylim([0 yl])
        line([x1 x1], [min(aspks_epoch(:,3)) yl], 'color','r','linestyle','--')
        line([x2 x2], [min(aspks_epoch(:,3)) yl], 'color','r','linestyle','--')
        axis on
        for s=1:size(aspks_epoch,1)
            if ~clusters_A.isMF(aspks_epoch(s,2)), cl = hex2rgb('d86d35'); else, cl = hex2rgb('7e4b8e'); end
            x = (aspks_epoch(s,1) - thisRip.STtime)*1e3+mar/2;
            patch([x-ti x+ti x+ti x-ti], [aspks_epoch(s,end)+.2 aspks_epoch(s,end)+.2 aspks_epoch(s,end)+.8 aspks_epoch(s,end)+.8],...
                cl,'edgecolor',cl,'edgealpha',0)
        end
        yticks([1.5:yl+.5])
        yticklabels(c.ID)
        axis ij

        xlabel(num2str(thisRip.RippleDuration*1000))
               %% save fig


        ROOT.Fig_en = [ROOT.Fig '\' suf1];
        if ~exist(ROOT.Fig_en), mkdir(ROOT.Fig_en); end
        saveas(gca,[ROOT.Fig_en '\'  thisRip.ID{1} '.svg'])
        saveas(gca,[ROOT.Fig_en '\'  thisRip.ID{1} '.png'])
        close all
catch
    close all
end

end
end

%%
function  [spks_epoch,spks_epoch_in,Units,UnitsA] = LoadUnits(thisRip,clusters,Spike,Params,mar)
spks_epoch=[]; spks_epoch_in=[];
u=0; Units={};UnitsA={};
if ismember('onMazeMaxFR_field', clusters.Properties.VariableNames)
    col4 = 'onMazeMaxFR_field';
else
    col4 = 'onMazeMaxFR';
end

cls_all = size(clusters,1);
for un = 1:cls_all
    [thisRID,thisSID,thisTTID,thisCLID, thisFLID] = parsing_clusterID(clusters.ID{un},1);

    Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);

    thisSpks = Spk.t_spk(Spk.t_spk>=thisRip.STtime-mar/Params.Fs & Spk.t_spk<=thisRip.EDtime+mar/Params.Fs);
    thisSpks_in = Spk.t_spk(Spk.t_spk>=thisRip.STtime & Spk.t_spk<=thisRip.EDtime);
    if ~isempty(thisSpks_in)

        s1 = ones(size(thisSpks,1),1); s2 = ones(size(thisSpks_in,1),1);
        %                 spks_epoch = [spks_epoch;[thisSpks,s1*un,s1*u,s1*clusters_B.SI(un),...
        %                     s1*clusters_B.RDI_LScene(un), s1*clusters_B.RDI_RScene(un), s1*clusters_B.RDI_LR(un), s1*clusters_A.(['RDI_hetero_' filter_ns])(h)]];
        spks_epoch = [spks_epoch;[thisSpks,s1*un,s1*u,s1*clusters.(col4)(un),...
            s1*clusters.RDI_LScene(un), s1*clusters.RDI_RScene(un), s1*clusters.RDI_LR(un),...
            s1*clusters.Selectivity_LScene(un),s1*clusters.Selectivity_RScene(un),s1*clusters.Selectivity_LR(un)]];
        spks_epoch_in = [spks_epoch_in;[thisSpks_in,s2*un,s2*u,s2*clusters.(col4)(un),...
            s2*clusters.RDI_LScene(un), s2*clusters.RDI_RScene(un), s2*clusters.RDI_LR(un),...
            s2*clusters.Selectivity_LScene(un),s2*clusters.Selectivity_RScene(un),s2*clusters.Selectivity_LR(un)]];
        % spks_epoch_in = [spks_epoch_in;[thisSpks_in,s2*un,s2*u,s2*clusters_B.SI(un),...
        %                     s2*clusters_B.RDI_LScene(un), s2*clusters_B.RDI_RScene(un), s2*clusters_B.RDI_LR(un),s2*clusters_A.(['RDI_hetero_' filter_ns])(h)]];
        u=u+1;
        if ~isnan(str2double(thisFLID))
            Units = [Units; [thisTTID '-' thisCLID '-' thisFLID]];
        else
            Units = [Units; [thisTTID '-' thisCLID]];
        end
        UnitsA = [UnitsA; [clusters.ID(un)]];
    end

end

% [~,ia,~] = unique(spks_epoch(:,3),'rows');
% spks_epoch_u = spks_epoch(ia,:);
% [~,ia,~] = unique(spks_epoch_in(:,3),'rows');
% spks_epoch_u_in = spks_epoch_in(ia,:);
end