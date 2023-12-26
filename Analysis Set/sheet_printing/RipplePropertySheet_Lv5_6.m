Initial_SWRFilter_common;
warning off


RegionList = {'CA1','SUB'};
for reg=1:2
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
ROOT.Fig = [ROOT.Processed '\ripples_mat\ProfilingSheet\R31_' thisRegion];
ROOT.Units = [ROOT.Processed '\units_mat\U2'];
ROOT.Behav = [ROOT.Processed '\behavior_mat'];

dir = '';

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


%%

RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis.xlsx']);
% UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion2 '.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_B;
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


%%

%%
% Isin=[];
% for uid = 1:size(UnitsTable_A,1)
%     id = UnitsTable_A.ID{uid};
%     for u=1:size(UnitsTable_B,1)
%         id2 = UnitsTable_B.ID{u};
%         k= strncmp(id,id2,12);
%         if k, break; end
%     end
%     Isin(uid) = k;
% end
% UnitsTable_A(~Isin,:)=[];
%
% writetable(UnitsTable_A,[ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx'])
% writetable(UnitsTable_B,[ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx'])

%%
% RipplesTable = sortrows(RipplesTable,{'nFields','Decoding_Rsq'},{'descend','ascend'});
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

            FRMaps_B = LoadFRMap(ROOT,clusters_B);
            FRMaps_A = LoadFRMap(ROOT,clusters_A);
            stem_end_index = stem_end_index*size(FRMaps_B,2)/48;
        end

        Clist=jet(256);
        %% Load units
        [spks_epoch,spks_epoch_in,spks_epoch_u,spks_epoch_u_in,Units,UnitsA] = LoadUnits(thisRip,clusters_B,Spike,Params,mar);

        [aspks_epoch,aspks_epoch_in,aspks_epoch_u,aspks_epoch_u_in,Unitsa,UnitsAa] = LoadUnits(thisRip,clusters_A,Spike,Params,mar);

        %% Set order
        [spks_epoch,spks_epoch_in,spks_epoch_u,spks_epoch_u_in,Units,UnitsA,FRMap_B,ord,flag] = ...
            SetOrder(spks_epoch,spks_epoch_in,spks_epoch_u,spks_epoch_u_in,Units,UnitsA,FRMaps_B);

        [aspks_epoch,aspks_epoch_in,aspks_epoch_u,aspks_epoch_u_in,Unitsa,UnitsAa,FRMap_A,orda,flag] = ...
            SetOrder(aspks_epoch,aspks_epoch_in,aspks_epoch_u,aspks_epoch_u_in,Unitsa,UnitsAa,FRMaps_A);
        %%
        %         size(aspks_epoch_u,1)
        if ~flag, continue; end

        %% smooth FRMap
        [FRMapsm_A,FRMapsm_L,FRMapsm_R,FRMapsm_C,start_index,end_index] = SmoothFRMap(FRMap_B);

        [FRMapsm_Aa,~,~,~,~,~] = SmoothFRMap(FRMap_A);

        %% load replay
        Replay = load([ROOT.Rip4 '\' RipplesTable.ID{sid} '.mat']);
        thisRip.Decoding_Rsq = Replay.R_actual(1); RipplesTable.Decoding_Rsq(sid)= Replay.R_actual(1);
        thisRip.DecodingDir = Replay.v{1}; RipplesTable.DecodingDir(sid) = Replay.v{1};
        %%
        figure('position',[317,63,1400,1000],'color','w');
        col=5;
        % title
        subplot(9,6,2:3)
        title(['(' thisRegion0 ') ' cell2mat(thisRip.ID) ', ' Exper ', ' dir],...
            'fontsize',15,'Interpreter','none')
        axis off

        % trial info
        subplot(9,6,4:5)
        if strcmp(Exper,'LSM'), CxtList = {'Zebra','Pebbles','Bamboo','Mountains'};
        elseif strcmp(Exper, 'SEB'), CxtList = {'Dot','Square','Zebra','Pebbles'};
        elseif strcmp(Exper, 'JS'), CxtList = {'Forest','','City'};
            if thisRip.area==5, thisRip.area=0; end
        end
        cxt = CxtList{thisRip.context};
        if thisRip.correctness==1, corr='Correct'; else, corr='Wrong'; end
        title(['trial ' (thisRip.trial{1}(end-2:end)) ', ' cxt ', ' corr],'fontsize',15)
        axis off

        % trial info
        subplot(9,6,6:6)
        %              title((['MultiVar cell ' jjnum2str(thisRip.(['nRDI_hetero_' filter_ns])*100,2) '%, mean RHS ' jjnum2str(thisRip.(['mRDI_hetero_' filter_ns]),2)]),'fontsize',15)
        xlim([1, size(spks_epoch_u_in,1)*0.7]);
        [~,b] = unique(spks_epoch_u_in(:,1));
        het = spks_epoch_u(b,[1,8]);
        for h=1:size(het,1)
            %                  text(h*1.2,0.5,jjnum2str(het(h,2),2))
        end
        axis off
        %% EEG
        thisEEG = EEG.(['TT' num2str(TargetTT_p)]).Raw(Ist:Ied);
        thisEEG_Ripple = EEG.(['TT' num2str(TargetTT_p)]).Filtered(Ist:Ied);
        hold on

        subplot(9,col,[1 col+1])
        plot(thisEEG,'k')
        title(['TT' num2str(TargetTT_p)])
        x1=mar; x2=mar+thisRip.RippleDuration*Params.Fs;
        dur2=x2+mar;
        xlim([0 dur2])
        line([x1 x1], [min(thisEEG) max(thisEEG)+50], 'color','r','linestyle','--')
        line([x2 x2], [min(thisEEG) max(thisEEG)+50], 'color','r','linestyle','--')
        axis off

        subplot(9,col,[col*2+1])
        plot(thisEEG_Ripple,'b')

        x1=mar; x2=mar+thisRip.RippleDuration*Params.Fs;
        dur2=x2+mar;
        xlim([0 dur2])
        line([x1 x1], [min(thisEEG_Ripple) max(thisEEG_Ripple)+50], 'color','r','linestyle','--')
        line([x2 x2], [min(thisEEG_Ripple) max(thisEEG_Ripple)+50], 'color','r','linestyle','--')
        axis off

        %% color scheme
        s = size(spks_epoch,1);
        clunits = hsv(max(spks_epoch(s,end))+1);
        clunits = [ones(s,1),zeros(s,1) , zeros(s,1)];

        %         aclunits = hsv(max(aspks_epoch(s,end))+1);
        clunits_L = clunits; clunits_R = clunits; clunits_C = clunits;

        for s=1:size(spks_epoch_u,1)

            if spks_epoch_u(s,8)<1 | abs(spks_epoch_u(s,5))<0.1
                clunits_L(spks_epoch_u(s,end)+1,:) = 1;
            end
            %             if(spks_epoch_u(s,8)>3)
            %                 clunits_L(spks_epoch_u(s,end)+1,:) = 0.7;
            %             end

            if spks_epoch_u(s,9)<1 | abs(spks_epoch_u(s,6))<0.1
                clunits_R(spks_epoch_u(s,end)+1,:) = 1;
            end
            %             if(spks_epoch_u(s,9)>3)
            %                 clunits_R(spks_epoch_u(s,end)+1,:) = 0.7;
            %             end

            if spks_epoch_u(s,10)<1 | abs(spks_epoch_u(s,7))<0.1
                clunits_C(spks_epoch_u(s,end)+1,:) = 1;
            end
            %             if(spks_epoch_u(s,10)>3)
            %                 clunits_C(spks_epoch_u(s,end)+1,:) = 0.7;
            %             end
        end

        %% reactivation rasterplot
        subplot(9,col,[3 5]*col+1)
        hold on
        x1=mar/2; x2=thisRip.RippleDuration*1e3+mar/2;
        xlim([0 dur2/2])
        ylim([0 max(aspks_epoch(:,3))+1])
        line([x1 x1], [min(aspks_epoch(:,3)) max(aspks_epoch(:,3))+1], 'color','r','linestyle','--')
        line([x2 x2], [min(aspks_epoch(:,3)) max(aspks_epoch(:,3))+1], 'color','r','linestyle','--')
        axis off
        for s=1:size(aspks_epoch,1)
            x = (aspks_epoch(s,1) - thisRip.STtime)*1e3+mar/2;
            patch([x-ti x+ti x+ti x-ti], [aspks_epoch(s,end)+.2 aspks_epoch(s,end)+.2 aspks_epoch(s,end)+.8 aspks_epoch(s,end)+.8],...
                'k','edgecolor','k','edgealpha',0)
        end
        %% Bayesian decoding plotting
        if 1
            marR = mar/Params.Fs/Params.tbinDuration;
            ax1 = subplot(9,col,[6 8]*col+1);
            hold on; axis on
            [pbinN, tbinN] = size(Replay.posterior{1});
            imagesc(ax1,Replay.posterior{1});
            colormap(ax1,flipud(gray))

            plot([-marR tbinN+marR]+0.5,[stem_end_index stem_end_index]+0.5, 'k:');
            line([.5 .5], [0 pbinN], 'color','r','linestyle','--')
            line([tbinN tbinN], [0 pbinN], 'color','r','linestyle','--')

            xlim([-.25 tbinN+.5]);
            ylim([0 pbinN]+0.5);
            set(gca, 'xtick', [0 tbinN]+0.5, 'xticklabel', {'0' [jjnum2str(thisRip.RippleDuration*1e3,1) 'ms']}, 'Fontsize', 10);
            set(gca,  'ytick', [0 stem_end_index pbinN]+0.5, 'yticklabel', {'Stbox', 'Dv', 'Fw'}, 'Fontsize', 10);
            xlabel('time(ms)');
            ylabel('Decoded position');
            %             title(['Bayesian decoding']);

            % plot the linear fitting line which has the best goodness of fit
            c_temp = Replay.c{1}([find(Replay.v{1}>0, 1) find(Replay.v{1}<0,1)]);
            v_temp = Replay.v{1}([find(Replay.v{1}>0, 1) find(Replay.v{1}<0,1)]);

            x = [0 tbinN]+0.5;
            for i = 1 : length(c_temp)
                plot(x,v_temp(i) * x + c_temp(i), '-', 'color', 'k', 'linewidth', 1);
                temp = v_temp(i) * x + c_temp(i);

            end
            if temp(1)<0, temp(1)=0; end
            if temp(2)<0, temp(2)=0; end
            slope = abs(temp(1)-temp(2));

                text(mean(x),mean(v_temp(i) * x + c_temp(i))*1.1,['R^2=' jjnum2str(thisRip.Decoding_Rsq,3)], 'color', 'r')
            hold off;
        end
        %% Reactivated cells
        temp = sortrows(spks_epoch_u,11);
        for c=1:size(Units,1)
            if isnan(temp(c,8))
                Units{c,2}='k';
            else
                Units{c,2} = 'k';
            end
        end

        cmap0 = flipud(gray);
        cmap='jet';

        ax1 = subplot(9,col,[1 2]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMapsm_Aa(1,:,:))'),Unitsa,[0 stem_end_index size(FRMaps_A,2)],{' ', ' ',' '}, 1);
        line([stem_end_index stem_end_index],[0 size(FRMaps_A,3)+1],'color','w','linestyle','--')
        title('Ensemble for Spatial React.')

        ax2 = subplot(9,col,[1 2]*col+3);
        fig = popul_FRMap(flip(squeeze(FRMapsm_A(1,:,:))'),Units,[0 stem_end_index size(FRMaps_B,2)],{' ', ' ',' '}, 1);
        line([stem_end_index stem_end_index],[0 size(FRMaps_B,3)+1],'color','w','linestyle','--')
        title('Ensemble for Non-Spatial React.')
        for i=1:size(FRMapsm_A,3)

            text(52,i,jjnum2str(temp(i,4),2),'color','k')

        end

        %%
        FRMap_Norm=[];
        id = clunits_L(:,1)-clunits_L(:,2);
        FRMap = FRMapsm_L;
%         FRMap(:,:,find(~id)) = nan;

        ULabel_1 = {}; ULabel_2 = {};
        for i=1:size(FRMap,3)
             FRMap_Norm(:,:,i) = FRMap(:,:,i) / nanmax(nanmax(nanmax(FRMap(:,:,i))));
            if isnan(nanmax(nanmax(nanmax(FRMap(:,:,i)))))
                ULabel_1{i,1} = ''; ULabel_2{i,1} = '';
            else
            ULabel_1{i,1} = jjnum2str(nanmax(nanmax(nanmax(FRMap(1,:,i)))),2);
            ULabel_2{i,1} = jjnum2str(nanmax(nanmax(nanmax(FRMap(2,:,i)))),2);
            end
        end


        ax3 = subplot(9,col,[3 4]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMap_Norm(1,:,:))'),ULabel_1 ,[0 stem_end_index size(FRMaps_B,2)],{' ', ' ',' '}, 1);
        line([stem_end_index stem_end_index],[0 size(FRMaps_B,3)+1],'color','w','linestyle','--')
        title('Zebra')

        ax4 = subplot(9,col,[3 4]*col+3);
        fig = popul_FRMap(flip(squeeze(FRMap_Norm(2,:,:))'),ULabel_2 ,[0 stem_end_index size(FRMaps_B,2)],{' ', ' ',' '}, 1);

        line([stem_end_index stem_end_index],[0 size(FRMaps_B,3)+1],'color','w','linestyle','--')
        title('Bamboo')
        %%
        id = clunits_R(:,1)-clunits_R(:,2);
        FRMap = FRMapsm_R;
%         FRMap(:,:,find(~id)) = nan;

   
        ULabel_1 = {}; ULabel_2 = {};
        for i=1:size(FRMap,3)
             FRMap_Norm(:,:,i) = FRMap(:,:,i) / nanmax(nanmax(nanmax(FRMap(:,:,i))));
            if isnan(nanmax(nanmax(nanmax(FRMap(:,:,i)))))
                ULabel_1{i,1} = ''; ULabel_2{i,1} = '';
            else
            ULabel_1{i,1} = jjnum2str(nanmax(nanmax(nanmax(FRMap(1,:,i)))),2);
            ULabel_2{i,1} = jjnum2str(nanmax(nanmax(nanmax(FRMap(2,:,i)))),2);
            end
        end

        ax5 = subplot(9,col,[5 6]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMap_Norm(1,:,:))'),ULabel_1,[0 stem_end_index size(FRMaps_B,2)],{' ', ' ',' '}, 1);
        line([stem_end_index stem_end_index],[0 size(FRMaps_B,3)+1],'color','w','linestyle','--')
        title('Pebbles')

        ax6 = subplot(9,col,[5 6]*col+3);
        fig = popul_FRMap(flip(squeeze(FRMap_Norm(2,:,:))'),ULabel_2,[0 stem_end_index size(FRMaps_B,2)],{' ', ' ',' '}, 1);
        line([stem_end_index stem_end_index],[0 size(FRMaps_B,3)+1],'color','w','linestyle','--')
        title('Mountain')

        %%
        id = clunits_C(:,1)-clunits_C(:,2);
        FRMap = FRMapsm_C;
%         FRMap(:,:,find(~id)) = nan;
        

    
        ULabel_1 = {}; ULabel_2 = {};
        for i=1:size(FRMap,3)
             FRMap_Norm(:,:,i) = FRMap(:,:,i) / nanmax(nanmax(nanmax(FRMap(:,:,i))));
            if isnan(nanmax(nanmax(nanmax(FRMap(:,:,i)))))
                ULabel_1{i,1} = ''; ULabel_2{i,1} = '';
            else
            ULabel_1{i,1} = jjnum2str(nanmax(nanmax(nanmax(FRMap(1,:,i)))),2);
            ULabel_2{i,1} = jjnum2str(nanmax(nanmax(nanmax(FRMap(2,:,i)))),2);
            end
        end

        ax7 = subplot(9,col,[7 8]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMap_Norm(1,1:stem_end_index,:))'),ULabel_1,[0 stem_end_index],{'Stbox', 'Dv'}, 1);
        line([stem_end_index stem_end_index],[0 size(FRMaps_B,3)+1],'color','w','linestyle','--')
        title('Left')

        ax8 = subplot(9,col,[7 8]*col+3);
        fig = popul_FRMap(flip(squeeze(FRMap_Norm(2,1:stem_end_index,:))'),ULabel_2,[0 stem_end_index],{'Stbox', 'Dv'}, 1);
        line([stem_end_index stem_end_index],[0 size(FRMaps_B,3)+1],'color','w','linestyle','--')
        title('Right')

        colormap(ax1,cmap)
        colormap(ax2,cmap)
        colormap(ax3,cmap)
        colormap(ax4,cmap)
        colormap(ax5,cmap)

        colormap(ax6,cmap)
        colormap(ax7,cmap)

        colormap(ax8,cmap)

        %% mean FR scatter
        % thisRip.pRDI_L=1; thisRip.pRDI_R=1; thisRip.pRDI_C=1;
        %         subplot(9,col,[3 4]*col+4);
        %         scatters_fr(FRMapsm_L, {CxtList{1},CxtList{3}},start_index,end_index,clunits_L,thisRip.pRDI_L_UV)
        %
        %         subplot(9,col,[5 6]*col+4);
        %         scatters_fr(FRMapsm_R, {CxtList{2},CxtList{4}},start_index,end_index,clunits_R,thisRip.pRDI_R_UV)
        %
        %         subplot(9,col,[7 8]*col+4);
        %         scatters_fr(FRMapsm_C, {'Left','Right'},start_index,end_index,clunits_C,thisRip.pRDI_C_UV)

        %% RDI bar graph
        %         subplot(9,col,[3 4]*col+5);
        subplot(9,col,[3 4]*col+4);
        bar_RDI(spks_epoch_u(:,end)+1,spks_epoch_u(:,5),clunits_L)
        if thisRip.nRDI_L_max < 5, thisRip.pRDI_L_max=nan; end
        title(['mean=' jjnum2str(thisRip.mRDI_L_max,3)])
        set ( gca, 'xdir', 'reverse' )


        %         subplot(9,col,[5 6]*col+5);
        subplot(9,col,[5 6]*col+4);
        bar_RDI(spks_epoch_u(:,end)+1,spks_epoch_u(:,6),clunits_R)
if thisRip.nRDI_R_max < 5, thisRip.pRDI_R_max=nan; end
           title(['mean=' jjnum2str(thisRip.mRDI_R_max,3)])
     set ( gca, 'xdir', 'reverse' )


        %         subplot(9,col,[7 8]*col+5);
        subplot(9,col,[7 8]*col+4);
        bar_RDI(spks_epoch_u(:,end)+1,spks_epoch_u(:,7),clunits_C)
if thisRip.nRDI_C_max < 5, thisRip.pRDI_C_max=nan; end
          title(['mean=' jjnum2str(thisRip.mRDI_C_max,3)])
       set ( gca, 'xdir', 'reverse' )

        %% permutation distribution
try
            S_perm = load([ROOT.Rip4 '\' thisRip.ID{1} '.mat'],"c","v",'R_shuffled','R_actual','posterior','p_test');

        subplot(9,col,[1 2]*col+5);
        DistPerm(S_perm.R_shuffled{1},S_perm.R_actual(1),thisRip.DecodingP_all,[0 .5],'r^2',0.5)
        title('permutation test')
        ylabel('r^2 for pos decoding')
end
         
   try
       NS_perm = load([ROOT.Rip5 ['\' thisRegion0 '-' thisRip.ID{1} '_max.mat']]);
        subplot(9,col,[3 4]*col+5);
        DistPerm(NS_perm.RDI_L.dist.mean,thisRip.mRDI_L_max,NS_perm.RDI_L.p_mean,[-1 1],'mean',0.2)
        ylabel('left scene RDI mean')
   end
   try
        subplot(9,col,[5 6]*col+5);
        DistPerm(NS_perm.RDI_R.dist.mean,thisRip.mRDI_R_max,NS_perm.RDI_R.p_mean,[-1 1],'mean',0.2)
        ylabel('right scene RDI mean')
   end
   try
        subplot(9,col,[7 8]*col+5);
        DistPerm(NS_perm.RDI_C.dist.mean,thisRip.mRDI_C_max,NS_perm.RDI_C.p_mean,[-1 1],'mean',0.2)
        ylabel('choice RDI mean')
end
        %% save fig
        if thisRip.DecodingP_all<0.05, suf1='Replay'; else, suf1='x'; end

%         if nanmin([thisRip.pRDI_L_UV, thisRip.pRDI_R_UV,thisRip.pRDI_C_UV])<0.05, suf2=['NSselective'];
if nanmin([thisRip.pRDI_L_max, thisRip.pRDI_R_max,thisRip.pRDI_C_max])<0.05
    suf2=['NSselective'];
else
    suf2='x';
end

        ROOT.Fig_en = [ROOT.Fig '\' suf1 '_' suf2];
        if ~exist(ROOT.Fig_en), mkdir(ROOT.Fig_en); end
        saveas(gca,[ROOT.Fig_en '\'  thisRip.ID{1} '.svg'])
        saveas(gca,[ROOT.Fig_en '\'  thisRip.ID{1} '.png'])
        close all

    catch
        close all
    end
end
writetable(RipplesTable,[ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis.xlsx'],'writemode','replacefile')
end

%%
function scatters_fr(FRMap, labels,start_index,end_index,clunits,p)
x=[]; y=[];
for i=1:size(FRMap,3)
    if start_index(i)
        x(i)= squeeze(nanmean(FRMap(1,start_index(i):end_index(i),i),2));
        y(i)= squeeze(nanmean(FRMap(2,start_index(i):end_index(i),i),2));
    else
        x(i)= 0;
        y(i)=0;
    end
end
l=length(x);
temp = [x',ones(l,1)+rand(l,1)*0.5-0.25 ; y',ones(l,1) *2.5+rand(l,1)*0.5-0.25];
c = [1:l,1:l]';
scatter(temp(:,2),temp(:,1),40,clunits(c,:),'filled','markerfacealpha',0.8)
hold on
for i=1:l
    line([temp(i,2),temp(i+l,2)],[temp(i,1),temp(i+l,1)],'color',clunits(i,:),'linewidth',.2)
end
xlim([0 3]); xticks([1 2.5]); xticklabels(labels)
ylim([0.01 1])
ylabel('mean field FR (Norm.)')
set(gca,'fontsize',8,'fontweight','b')
% [h,p] = ttest2(x,y);
%  if p<0.005
%      s='***';
%  elseif p<0.01
%      s='**';
%  elseif p<0.05
%      s='*';
%  else
%      s='n.s.';
%  end

if p<0.001
    s='<0.001';
else
    s=jjnum2str(p,3);
end

text(1.75,.8,s,'fontsize',10,'fontweight','b')
end
%%
function bar_RDI(x,y,clunits)
b= barh(x,y);
b.FaceColor = 'flat';
for c=1:size(b.CData,1)
    b.CData(c,:)= clunits(c,:);
end
xlim([-1.2 1.2]);
yticks([1:max(x)+1]); yticklabels({})
set(gca,'fontsize',8,'fontweight','b')
line([-.1 -.1],[0 max(x)+1],'color','r','linestyle',':')
line([.1 .1],[0 max(x)+1],'color','r','linestyle',':')
end


%%
function  [spks_epoch,spks_epoch_in,spks_epoch_u,spks_epoch_u_in,Units,UnitsA] = LoadUnits(thisRip,clusters,Spike,Params,mar)
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

[~,ia,~] = unique(spks_epoch(:,3),'rows');
spks_epoch_u = spks_epoch(ia,:);
[~,ia,~] = unique(spks_epoch_in(:,3),'rows');
spks_epoch_u_in = spks_epoch_in(ia,:);
end

function [spks_epoch,spks_epoch_in,spks_epoch_u,spks_epoch_u_in,Units,UnitsA,FRMap,ord,flag] = SetOrder(spks_epoch,spks_epoch_in,spks_epoch_u,spks_epoch_u_in,Units,UnitsA,FRMaps_B)
flag=1;
if ~isempty(UnitsA)
    o0=[]; o1=[]; o2=[];
    o0=[1:length(UnitsA)]';
    [~,ord]=sort(spks_epoch_u_in(o0,1),'descend');
    UnitsC = UnitsA(ord);
else
    flag=0;
end
if isempty(ord), flag=0; end
if length(ord)<=1, flag=0; end
if flag
    %%

    for s=1:size(spks_epoch,1)
        spks_epoch(s,11) = find(ord==spks_epoch(s,3)+1)-1;
    end
    [~,ia] = sort(spks_epoch(:,11));
    spks_epoch = spks_epoch(ia,:);
    [~,ia,~] = unique(spks_epoch(:,3),'rows');
    spks_epoch_u = spks_epoch(ia,:);
    Units = Units(ord);
    UnitsA=UnitsA(ord);

    id = spks_epoch_u(:,2);
    FRMap=FRMaps_B(:,:,id);
    FRMap=FRMap(:,:,ord);
    if ~max(max(max(~isnan(FRMap(:,43:end,:)))))
        FRMap=FRMap(:,1:42,:);
    end
end
end
%%
function [FRMapsm_A,FRMapsm_L,FRMapsm_R,FRMapsm_C,start_index,end_index] = SmoothFRMap(FRMap)
FRMap_A = FRMap(1,:,:);FRMap_L=FRMap(2:3,:,:);FRMap_R=FRMap(4:5,:,:);FRMap_C=FRMap(6:7,:,:);
FRMapsm_A=[]; FRMapsm_L=[]; FRMapsm_R=[]; FRMapsm_C=[];
start_index=[]; end_index = [];
for i=1:size(FRMap,3)
    try
        [field_count, start_index(i), end_index(i), field_size] = field_boundary_function_jm(FRMap(1,:,i),0.25);
    catch
        start_index(i) = 0;
        end_index(i)=0;
    end
                FRMap_A(:,:,i)  = FRMap_A(:,:,i)/max(max(max(FRMap_A(:,:,i))));
    %             FRMap_L(:,:,i)  = FRMap_L(:,:,i)/max(max(max(FRMap_L(:,:,i))));
    %             FRMap_R(:,:,i)  = FRMap_R(:,:,i)/max(max(max(FRMap_R(:,:,i))));
    %             FRMap_C(:,:,i)  = FRMap_C(:,:,i)/max(max(max(FRMap_C(:,:,i))));
    for j=1:2
        FRMapsm_L(j,:,i) = smooth(FRMap_L(j,:,i), "moving");
        FRMapsm_R(j,:,i) = smooth(FRMap_R(j,:,i), "moving");
        FRMapsm_C(j,:,i) = smooth(FRMap_C(j,:,i), "moving");
    end
    FRMapsm_A(1,:,i) = smooth(FRMap_A(1,:,i), "moving");
end
end
%%
function DistPerm(dist,act,p,lim,ylab,y1)
binw=(lim(2)-lim(1))/40;
hold on
histogram(dist,'binwidth',binw,'facecolor','k','Normalization','probability')
line([act act],[0 y1],'color','r')
line([mean(lim) mean(lim)],[0 y1],'color','k','linestyle',':')
xlim(lim)
if p<0.001
    text(0.3,y1*0.8,['perm.p<0.001'])
else
    text(0.3,y1*0.8,['perm.p= ' jjnum2str(p,3)])
end
text(0.3,y1*0.9,[ylab '= ' jjnum2str(act,3)],'color','r')
set(gca,'xdir','reverse','ydir','normal')
end
