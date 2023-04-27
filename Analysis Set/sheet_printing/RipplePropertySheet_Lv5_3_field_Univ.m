Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Mother '\Processed Data\ripples_mat\R4'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R11_sub_field_uv'];
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
plot_matrix
filter_ns = 'C';


%%
col_name_A = string(UnitsTable_A.Properties.VariableNames(24:32));
UnitsTable_B(1,28:36) = array2table(nan(1,9));
UnitsTable_B = renamevars(UnitsTable_B, UnitsTable_B.Properties.VariableNames(28:36),col_name_A);
for uid=1:size(UnitsTable_B,1)
    tar = find(strncmp(UnitsTable_A.ID,UnitsTable_B.ID(uid),12));
    UnitsTable_B(uid,28:36) = UnitsTable_A(tar,24:32);
end
%%
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

            %% remove heterogeneous cells

            clusters_B = clusters_B(~clusters_B.(['MultiVar_' filter_ns]),:);
            clusters_A = clusters_A(~clusters_A.(['MultiVar_' filter_ns]),:);
            %% Load FRMap

            FRMaps = [];
            for cl=1:size(clusters_B.ID,1)
                try
                    thisField = table;
                    temp=[];

                    thisTTID = num2str(str2double(clusters_B.ID{cl}(8:9)));
                    thisCLID =  clusters_B.ID{cl}(11:12);
                    thisFLID =  jmnum2str(str2double(clusters_B.ID{cl}(14:end)),2);

                    thisMap.thisFieldMap1D = load([ROOT.Raw.Map '\rat' thisRSID '-' thisTTID '-' thisCLID '-' thisFLID '.mat']);

                    for i=1:numel(thisMap.thisFieldMap1D.skaggsMap1D)
                        FRMaps(i,1:length(thisMap.thisFieldMap1D.skaggsMap1D{i}),cl) = thisMap.thisFieldMap1D.skaggsMap1D{i};
                    end

                    FRMaps(6,:,cl) = thisMap.thisFieldMap1D.skaggsMap_left1D;
                    FRMaps(7,:,cl) = thisMap.thisFieldMap1D.skaggsMap_right1D;

                    %             FRMaps(:,:,cl) = FRMaps(:,:,cl) / max(max(FRMaps(:,:,cl)));
                catch
                end


            end
            for cl=1:size(clusters_B.ID,1)
                if sum(FRMaps(1:7,:,cl),'all','omitnan')==0, FRMaps(1:7,:,cl)=nan; end
            end
            stem_end_index = stem_end_index*size(FRMaps,2)/48;
        end


        %% Load units
        spks_epoch=[]; spks_epoch_in=[];
        u=0; Units={};UnitsA={};
        Clist=jet(256);

        cls_all = size(clusters_B,1);
        for un = 1:cls_all
           [thisRID,thisSID,thisTTID,thisCLID, thisFLID] = parsing_clusterID(clusters_B.ID{un},1);
            
            Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);

            thisSpks = Spk.t_spk(Spk.t_spk>=thisRip.STtime-mar/Params.Fs & Spk.t_spk<=thisRip.EDtime+mar/Params.Fs);
            thisSpks_in = Spk.t_spk(Spk.t_spk>=thisRip.STtime & Spk.t_spk<=thisRip.EDtime);
            if ~isempty(thisSpks_in)
                h = contains(clusters_A.ID,clusters_B.ID{un}(1:12));
                s1 = ones(size(thisSpks,1),1); s2 = ones(size(thisSpks_in,1),1);
                spks_epoch = [spks_epoch;[thisSpks,s1*un,s1*u,s1*clusters_B.SI(un),...
                    s1*clusters_B.RDI_LScene(un), s1*clusters_B.RDI_RScene(un), s1*clusters_B.RDI_LR(un), s1*clusters_A.(['RDI_hetero_' filter_ns])(h)]];
                spks_epoch_in = [spks_epoch_in;[thisSpks_in,s2*un,s2*u,s2*clusters_B.SI(un),...
                    s2*clusters_B.RDI_LScene(un), s2*clusters_B.RDI_RScene(un), s2*clusters_B.RDI_LR(un),s2*clusters_A.(['RDI_hetero_' filter_ns])(h)]];
                u=u+1;
                Units = [Units; [thisTTID '-' thisCLID '-' thisFLID]];
                UnitsA = [UnitsA; [clusters_B.ID(un)]];
            end

        end

        [~,ia,~] = unique(spks_epoch(:,3),'rows');
        spks_epoch_u = spks_epoch(ia,:);
        [~,ia,~] = unique(spks_epoch_in(:,3),'rows');
        spks_epoch_u_in = spks_epoch_in(ia,:);



        %% Set order
        if ~isempty(UnitsA)
            o0=[]; o1=[]; o2=[];
            o0=[1:length(UnitsA)]';
            [~,ord]=sort(spks_epoch_u_in(o0,1),'descend');
            UnitsC = UnitsA(ord);
        else
            continue;
        end
        if isempty(ord), continue; end
        if length(ord)<=1, continue; end
        %%

        for s=1:size(spks_epoch,1)
            spks_epoch(s,9) = find(ord==spks_epoch(s,3)+1)-1;
        end
        [~,ia] = sort(spks_epoch(:,9));
        spks_epoch = spks_epoch(ia,:);
        [~,ia,~] = unique(spks_epoch(:,3),'rows');
        spks_epoch_u = spks_epoch(ia,:);
        Units = Units(ord);
        UnitsA=UnitsA(ord);

        id = spks_epoch_u(:,2);
        FRMap=FRMaps(:,:,id);
        FRMap=FRMap(:,:,ord);
        if ~max(max(max(~isnan(FRMap(:,43:end,:)))))
            FRMap=FRMap(:,1:42,:);
        end
        %% smooth FRMap
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
            FRMap_L(:,:,i)  = FRMap_L(:,:,i)/max(max(max(FRMap_L(:,:,i))));
            FRMap_R(:,:,i)  = FRMap_R(:,:,i)/max(max(max(FRMap_R(:,:,i))));
            FRMap_C(:,:,i)  = FRMap_C(:,:,i)/max(max(max(FRMap_C(:,:,i))));
            for j=1:2
                FRMapsm_L(j,:,i) = smooth(FRMap_L(j,:,i), "moving");
                FRMapsm_R(j,:,i) = smooth(FRMap_R(j,:,i), "moving");
                FRMapsm_C(j,:,i) = smooth(FRMap_C(j,:,i), "moving");
            end
            FRMapsm_A(1,:,i) = smooth(FRMap_A(1,:,i), "moving");
        end
        %% load replay
        Replay = load([ROOT.Rip4 '\' RipplesTable.ID{sid} '.mat']);
        %%
        figure('position',[317,63,1400,915],'color','w');
        col=5;
        % title
        subplot(9,6,2)
        title([cell2mat(thisRip.ID) ', ' Exper ', ' dir],'fontsize',15)
        axis off

        % trial info
        subplot(9,6,3:4)
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
        subplot(9,6,5:6)
             title((['MultiVar cell ' jjnum2str(thisRip.(['nRDI_hetero_' filter_ns])*100,2) '%, mean RHS ' jjnum2str(thisRip.(['mRDI_hetero_' filter_ns]),2)]),'fontsize',15)
             xlim([1, size(spks_epoch_u_in,1)*0.7]);
             [~,b] = unique(spks_epoch_u_in(:,1));
             het = spks_epoch_u(b,[1,8]);
             for h=1:size(het,1)
                 text(h*1.2,0.5,jjnum2str(het(h,2),2))
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

        %% reactivation rasterplot
        clunits = hsv(max(spks_epoch(s,end))+1);
        subplot(9,col,[3 5]*col+1)
        hold on
        x1=mar/2; x2=thisRip.RippleDuration*1e3+mar/2;
        xlim([0 dur2/2])
        ylim([0 max(spks_epoch(:,3))+1])
        line([x1 x1], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','r','linestyle','--')
        line([x2 x2], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','r','linestyle','--')
        axis off
        for s=1:size(spks_epoch,1)
            x = (spks_epoch(s,1) - thisRip.STtime)*1e3+mar/2;
            patch([x-ti x+ti x+ti x-ti], [spks_epoch(s,end)+.2 spks_epoch(s,end)+.2 spks_epoch(s,end)+.8 spks_epoch(s,end)+.8],...
                clunits(spks_epoch(s,end)+1,:),'edgecolor','k','edgealpha',0)
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
            if thisRip.DecodingP_all<0.001
                text(mean(x),mean(v_temp(i) * x + c_temp(i))*1.1,'<0.001')
            else
            text(mean(x),mean(v_temp(i) * x + c_temp(i))*1.1,jjnum2str(thisRip.DecodingP_all,2))
            end
            hold off;
        end
        %% Reactivated cells
        temp = sortrows(spks_epoch_u,9);
for c=1:size(Units,1)
    if isnan(temp(c,8))
        Units{c,2}='k';
    else
        Units{c,2} = 'k';
    end
end

        cmap = flipud(gray);
        cmap='jet';
        ax2 = subplot(9,col,[1 2]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMapsm_A(1,:,:))'),Units,[0 stem_end_index size(FRMap,2)],{' ', ' ',' '}, 1);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        title('Overall')

        ax3 = subplot(9,col,[3 4]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMapsm_L(1,:,:))'),Units,[0 stem_end_index size(FRMap,2)],{' ', ' ',' '}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        title('Zebra')

        ax4 = subplot(9,col,[3 4]*col+3);
        fig = popul_FRMap(flip(squeeze(FRMapsm_L(2,:,:))'),Units,[0 stem_end_index size(FRMap,2)],{' ', ' ',' '}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        title('Bamboo')

        ax5 = subplot(9,col,[5 6]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMapsm_R(1,:,:))'),Units,[0 stem_end_index size(FRMap,2)],{' ', ' ',' '}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        title('Pebbles')

        ax6 = subplot(9,col,[5 6]*col+3);
        fig = popul_FRMap(flip(squeeze(FRMapsm_R(2,:,:))'),Units,[0 stem_end_index size(FRMap,2)],{' ', ' ',' '}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        title('Mountain')

        ax7 = subplot(9,col,[7 8]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMapsm_C(1,1:stem_end_index,:))'),Units,[0 stem_end_index],{'Stbox', 'Dv'}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        title('Left')

        ax8 = subplot(9,col,[7 8]*col+3);
        fig = popul_FRMap(flip(squeeze(FRMapsm_C(2,1:stem_end_index,:))'),Units,[0 stem_end_index],{'Stbox', 'Dv'}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        title('Right')

        colormap(ax2,cmap)
        colormap(ax3,cmap)
        colormap(ax4,cmap)
        colormap(ax5,cmap)

        colormap(ax6,cmap)
        colormap(ax7,cmap)

        colormap(ax8,cmap)

        %% mean FR scatter
        % thisRip.pRDI_L=1; thisRip.pRDI_R=1; thisRip.pRDI_C=1;
        subplot(9,col,[3 4]*col+4);
        scatters_fr(FRMapsm_L, {CxtList{1},CxtList{3}},start_index,end_index,clunits,thisRip.pRDI_L_UV)

        subplot(9,col,[5 6]*col+4);
        scatters_fr(FRMapsm_R, {CxtList{2},CxtList{4}},start_index,end_index,clunits,thisRip.pRDI_R_UV)

        subplot(9,col,[7 8]*col+4);
        scatters_fr(FRMapsm_C, {'Left','Right'},start_index,end_index,clunits,thisRip.pRDI_C_UV)

        %% RDI bar graph
        subplot(9,col,[3 4]*col+5);
        bar_RDI(spks_epoch_u(:,end)+1,spks_epoch_u(:,5),clunits)
        set ( gca, 'xdir', 'reverse' )

        subplot(9,col,[5 6]*col+5);
        bar_RDI(spks_epoch_u(:,end)+1,spks_epoch_u(:,6),clunits)
        set ( gca, 'xdir', 'reverse' )

        subplot(9,col,[7 8]*col+5);
        bar_RDI(spks_epoch_u(:,end)+1,spks_epoch_u(:,7),clunits)
        set ( gca, 'xdir', 'reverse' )

        %% save fig
        if thisRip.DecodingP_all<0.05, suf1='Replay'; else, suf1='x'; end
        suf1='';
        if thisRip.(['pRDI_' filter_ns '_UV'])<0.05, suf2=[filter_ns '_selective']; else, suf2='x'; end

        ROOT.Fig_en = [ROOT.Fig '\' suf2];
        if ~exist(ROOT.Fig_en), mkdir(ROOT.Fig_en); end
        saveas(gca,[ROOT.Fig_en '\'  thisRip.ID{1} '.svg'])
        saveas(gca,[ROOT.Fig_en '\'  thisRip.ID{1} '.png'])
        close all

    catch
        close all
    end
end

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

if p<0.005
    s='***';
else
    s=jjnum2str(p,3);
end

text(1.75,.8,s,'fontsize',10,'fontweight','b')
end

function bar_RDI(x,y,clunits)
b= barh(x,y);
b.FaceColor = 'flat';
b.CData= clunits;
xlim([-1.2 1.2]);
yticks([1:max(x)+1]); yticklabels({})
set(gca,'fontsize',8,'fontweight','b')
end
