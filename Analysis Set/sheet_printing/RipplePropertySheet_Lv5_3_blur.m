Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data (Blur)'];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Save '\ripples_mat\R4'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R11_AI'];
ROOT.Units = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];

dir = '';
ROOT.Fig = [ROOT.Fig3];
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion '.xlsx']);
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
%%
for sid=1:size(RipplesTable,1)
    try
        thisRip = RipplesTable(sid,:);
        %         if thisRip.rat~=80, continue; end
        thisRID = jmnum2str(thisRip.rat,3);
        thisSID = jmnum2str(thisRip.session,2);
%         if strcmp(thisRID,'561'), continue; end
   


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
            TargetTT = GetTargetTT(ROOT,thisRSID,thisRegion,Params,0);

            thisTT_table = TT_table(TT_table.rat==thisRip.rat & TT_table.session==thisRip.session,:);
%             for t=1:size(thisTT_table,1)
%                 if ~ismember(thisTT_table.TT(t),TargetTT)
%                     thisTT_table.TT(t)=0;
%                 end
%             end
%             thisTT_table= thisTT_table(thisTT_table.TT~=0,:);
%             [~,t] = max(thisTT_table.RippleBandMean);
%             TargetTT_p = thisTT_table.TT(t);
TargetTT_p = TargetTT(1);

            EEG = LoadEEGData(ROOT, thisRSID, TargetTT_p,Params,Params_Ripple);

            %             Pos = load([ROOT.Raw.Mother '\rat' thisRSID(1:3) '\rat' thisRSID '\ParsedPosition.mat']);
            clusters = UnitsTable(UnitsTable.rat==thisRip.rat & UnitsTable.session==thisRip.session,:);
            clusters_A = UnitsTable_A(UnitsTable_A.rat==thisRip.rat & UnitsTable_A.session==thisRip.session,:);
            Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);
            load([ROOT.Rip0 '\' thisRSID '.mat'])
            disp([thisRSID ' plotting...'])
            thisRSID_old = thisRSID;

            %% Load FRMap

        FRMaps = [];
        for cl=1:size(clusters_A.ID,1)
            try
            thisField = table;
            temp=[];

            thisTTID = num2str(str2double(clusters_A.ID{cl}(8:9)));
            thisCLID =  clusters_A.ID{cl}(11:12);

            thisMap.thisFieldMap1D = load([ROOT.Raw.Map '\rat' thisRSID '-' thisTTID '-' thisCLID  '.mat']);

            for i=1:numel(thisMap.thisFieldMap1D.skaggsMap1D)
                FRMaps(i,1:length(thisMap.thisFieldMap1D.skaggsMap1D{i}),cl) = thisMap.thisFieldMap1D.skaggsMap1D{i};
            end

            FRMaps(11,:,cl) = thisMap.thisFieldMap1D.skaggsMap_left1D;
            FRMaps(12,:,cl) = thisMap.thisFieldMap1D.skaggsMap_right1D;

%             FRMaps(:,:,cl) = FRMaps(:,:,cl) / max(max(FRMaps(:,:,cl)));
            catch
            end


        end
        for cl=1:size(clusters_A.ID,1)
            if sum(FRMaps(1:12,:,cl),'all','omitnan')==0, FRMaps(1:12,:,cl)=nan; end
        end
        stem_end_index = stem_end_index*size(FRMaps,2)/48;
        end


        %% Load units
        spks_epoch=[]; spks_epoch_in=[];
        u=0; Units={};UnitsA={};
        Clist=jet(256);

        clusters_A =  clusters;
%         clusters_A = sortrows(clusters_A,thisRDI,'ascend');
%         unz = find(clusters_A.(thisRDI)<0,1,'last');

        cls_all = size(clusters_A,1);
        for un = 1:cls_all
            thisTTID = num2str(clusters_A.TT(un));
            thisCLID = num2str(str2double(clusters_A.ID{un}(end-1:end)));
            Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);

            thisSpks = Spk.t_spk(Spk.t_spk>=thisRip.STtime-mar/Params.Fs & Spk.t_spk<=thisRip.EDtime+mar/Params.Fs);
            thisSpks_in = Spk.t_spk(Spk.t_spk>=thisRip.STtime & Spk.t_spk<=thisRip.EDtime);
            if ~isempty(thisSpks_in)
                s1 = ones(size(thisSpks,1),1); s2 = ones(size(thisSpks_in,1),1);
                spks_epoch = [spks_epoch;[thisSpks,s1*un,s1*u,s1*clusters_A.SI(un),...
                    s1*clusters_A.RDI_No(un), s1*clusters_A.RDI_Lo(un), s1*clusters_A.RDI_Hi(un)]];
                spks_epoch_in = [spks_epoch_in;[thisSpks_in,s2*un,s2*u,s2*clusters_A.SI(un),...
                    s2*clusters_A.RDI_No(un), s2*clusters_A.RDI_Lo(un), s2*clusters_A.RDI_Hi(un)]];
                u=u+1;
                Units = [Units; [thisTTID '-' thisCLID]];
                UnitsA = [UnitsA; [clusters_A.ID(un)]];
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
            spks_epoch(s,8) = find(ord==spks_epoch(s,3)+1)-1;
        end
        [~,ia] = sort(spks_epoch(:,8));
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
       FRMap_A = FRMap(10,:,:);FRMap_No=FRMap([4 7],:,:);FRMap_Lo=FRMap([5 8],:,:);FRMap_Hi=FRMap([6 9],:,:);
       FRMapsm_A=[]; FRMapsm_No=[]; FRMapsm_Lo=[]; FRMapsm_Hi=[];
       start_index=[]; end_index = [];
        for i=1:size(FRMap,3)
            try
            [field_count, start_index(i), end_index(i), field_size, h] = field_boundary_function_jm(FRMap(1,:,i),0.25);
            catch
                start_index(i) = 0;
                end_index(i)=0;
            end
            FRMap_A(:,:,i)  = FRMap_A(:,:,i)/max(max(max(FRMap_A(:,:,i))));
            FRMap_No(:,:,i)  = FRMap_No(:,:,i)/max(max(max(FRMap_No(:,:,i))));
            FRMap_Lo(:,:,i)  = FRMap_Lo(:,:,i)/max(max(max(FRMap_Lo(:,:,i))));
            FRMap_Hi(:,:,i)  = FRMap_Hi(:,:,i)/max(max(max(FRMap_Hi(:,:,i))));
            for j=1:2
                FRMapsm_No(j,:,i) = smooth(FRMap_No(j,:,i), "moving");
                FRMapsm_Lo(j,:,i) = smooth(FRMap_Lo(j,:,i), "moving");
                FRMapsm_Hi(j,:,i) = smooth(FRMap_Hi(j,:,i), "moving");
            end
            FRMapsm_A(1,:,i) = smooth(FRMap_A(1,:,i), "moving");
        end
        %% load replay
        Replay = load([ROOT.Rip4 '\' RipplesTable.ID{sid} '.mat']);
        %%
        figure('position',[317,63,1400,915],'color','w');
        col=5;
        % title
        subplot(9,6,3)
        title([cell2mat(thisRip.ID) ', ' Exper ', ' dir],'fontsize',15)
        axis off

        % trial info
        subplot(9,6,6)
        if strcmp(Exper,'LSM'), CxtList = {'Zebra','Pebbles','Bamboo','Mountains'};
        elseif strcmp(Exper, 'SEB'), CxtList = {'Dot','Square','Zebra','Pebbles'};
        elseif strcmp(Exper, 'JS'), CxtList = {'Forest','','City'};
            if thisRip.area==5, thisRip.area=0; end
        end
        cxt = CxtList{thisRip.context};
        if thisRip.correctness==1, corr='Correct'; else, corr='Wrong'; end
        if thisRip.ambiguity==1
           if str2double(thisRip.trial{1}(end-2:end))<=20
               amb='STD';
           else, amb='No-Blur';
           end
        elseif thisRip.ambiguity==2
            amb='Low-Blur';
        elseif thisRip.ambiguity==3
            amb='High-Blur';
        else
            amb='Blur Unknown';
        end
        title(['trial ' (thisRip.trial{1}(end-2:end)) ', ' cxt ', ' amb ', ' corr],'fontsize',15)
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
        clunits = hsv(max(spks_epoch(s,8))+1);
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
            patch([x-ti x+ti x+ti x-ti], [spks_epoch(s,8)+.2 spks_epoch(s,8)+.2 spks_epoch(s,8)+.8 spks_epoch(s,8)+.8],...
                clunits(spks_epoch(s,8)+1,:),'edgecolor','k','edgealpha',0)
        end
        %% Bayesian decoding plotting
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
        text(mean(x),mean(v_temp(i) * x + c_temp(i))*1.1,jjnum2str(thisRip.DecodingP_all,2))


        hold off;

        %% Reactivated cells
        
        cmap = flipud(gray);
        cmap='jet';
        ax2 = subplot(9,col,[1 2]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMapsm_A(1,:,:))'),Units,[0 stem_end_index size(FRMap,2)],{' ', ' ',' '}, 1);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
         title('Overall')

        ax3 = subplot(9,col,[3 4]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMapsm_No(1,:,:))'),Units,[0 stem_end_index size(FRMap,2)],{' ', ' ',' '}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        title('No Blur-Zebra(L)')

        ax4 = subplot(9,col,[3 4]*col+3);
        fig = popul_FRMap(flip(squeeze(FRMapsm_No(2,:,:))'),Units,[0 stem_end_index size(FRMap,2)],{' ', ' ',' '}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
         title('No Blur-Pebble(R)')

        ax5 = subplot(9,col,[5 6]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMapsm_Lo(1,:,:))'),Units,[0 stem_end_index size(FRMap,2)],{' ', ' ',' '}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
         title('Low Blur-Zebra(L)')

        ax6 = subplot(9,col,[5 6]*col+3);
        fig = popul_FRMap(flip(squeeze(FRMapsm_Lo(2,:,:))'),Units,[0 stem_end_index size(FRMap,2)],{' ', ' ',' '}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
         title('Low Blur-Pebble(R)')

        ax7 = subplot(9,col,[7 8]*col+2);
        fig = popul_FRMap(flip(squeeze(FRMapsm_Hi(1,1:stem_end_index,:))'),Units,[0 stem_end_index],{'Stbox', 'Dv'}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
        title('High Blur-Zebra(L)')

        ax8 = subplot(9,col,[7 8]*col+3);
        fig = popul_FRMap(flip(squeeze(FRMapsm_Hi(2,1:stem_end_index,:))'),Units,[0 stem_end_index],{'Stbox', 'Dv'}, 0);
        line([stem_end_index stem_end_index],[0 size(FRMap,3)+1],'color','w','linestyle','--')
         title('High Blur-Pebble(R)')

        colormap(ax2,cmap)
        colormap(ax3,cmap)
                colormap(ax4,cmap)
        colormap(ax5,cmap)

                colormap(ax6,cmap)
        colormap(ax7,cmap)

                colormap(ax8,cmap)

%% mean FR scatter
subplot(9,col,[3 4]*col+4);
scatters_fr(FRMapsm_No, {CxtList{1},CxtList{3}},start_index,end_index,clunits,thisRip.pRDI_No)

subplot(9,col,[5 6]*col+4);
scatters_fr(FRMapsm_Lo, {CxtList{2},CxtList{4}},start_index,end_index,clunits,thisRip.pRDI_Lo)

subplot(9,col,[7 8]*col+4);
scatters_fr(FRMapsm_Hi, {'Left','Right'},start_index,end_index,clunits,thisRip.pRDI_Hi)
        %% RDI bar graph
subplot(9,col,[3 4]*col+5);
bar_RDI(spks_epoch_u(:,8)+1,spks_epoch_u(:,5),clunits)

subplot(9,col,[5 6]*col+5);
bar_RDI(spks_epoch_u(:,8)+1,spks_epoch_u(:,6),clunits)

subplot(9,col,[7 8]*col+5);
bar_RDI(spks_epoch_u(:,8)+1,spks_epoch_u(:,7),clunits)
        %% save fig
                if thisRip.DecodingP_all<0.05, suf1='Replay'; else, suf1='x'; end
        if nanmin([thisRip.pRDI_All,thisRip.pRDI_No,thisRip.pRDI_Lo,thisRip.pRDI_Hi])<0.05, suf2='Bias'; else, suf2='x'; end
        
        ROOT.Fig_en = [ROOT.Fig '\' suf1 '_' suf2];
        if ~exist(ROOT.Fig_en), mkdir(ROOT.Fig_en); end
        saveas(gca,[ROOT.Fig_en '\' num2str(RipplesTable.nPCs(sid)) '_' thisRip.ID{1} '.svg'])
        saveas(gca,[ROOT.Fig_en '\' num2str(RipplesTable.nPCs(sid)) '_' thisRip.ID{1} '.png'])
        close all


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
 if p<0.005
     s='***';
 elseif p<0.01
     s='**';
 elseif p<0.05
     s='*';
 else
     s='n.s.';
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