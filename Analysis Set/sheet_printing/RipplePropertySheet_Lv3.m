Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R3'];
ROOT.Fig3 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R3'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
if ~exist(ROOT.Fig3), mkdir(ROOT.Fig3); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 0.5;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.1*Params.Fs;
dur = 0.4*Params.Fs;
%%
for sid=1:size(RipplesTable,1)
    try
        thisRip = RipplesTable(sid,:);
%         if thisRip.rat~=80, continue; end
        Ist = thisRip.STindex-mar;
        Ied = thisRip.EDindex+mar;
        temp_r = decimalToBinaryVector(thisRip.TTs);
        RippleTT = flip(length((temp_r))+1 - find(temp_r))';
        Exper = cell2mat(thisRip.experimenter);
        if ~ismember(Exper,Experimenter), continue; end
        thisRSID = [jmnum2str(thisRip.rat,3) '-' jmnum2str(thisRip.session,2)];
        
        
        if ~strcmp(thisRSID, thisRSID_old)
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
            
            EEG = LoadEEGData(ROOT, thisRSID, TargetTT,Params,Params_Ripple);
            NumTT = length(fieldnames(EEG));
            Pos = load([ROOT.Raw.Mother '\rat' thisRSID(1:3) '\rat' thisRSID '\ParsedPosition.mat']);
            clusters = UnitsTable(UnitsTable.rat==thisRip.rat & UnitsTable.session==thisRip.session,:);
            clusters_A = UnitsTable_A(UnitsTable_A.rat==thisRip.rat & UnitsTable_A.session==thisRip.session,:);
            Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);
            load([ROOT.Rip0 '\' thisRSID '.mat'])
            disp([thisRSID ' plotting...'])
            thisRSID_old = thisRSID;
        end
        %%
        figure('position',[317,63,923,915],'color','w');
        % title
        subplot(9,6,1)
        title([cell2mat(thisRip.ID) ', ' Exper],'fontsize',15)
        axis off
        
        % trial info
        subplot(9,6,4)
        if strcmp(Exper,'LSM'), CxtList = {'Zebra','Pebbles','Bamboo','Mountains'};
        elseif strcmp(Exper, 'SEB'), CxtList = {'Dot','Square','Zebra','Pebbles'};
        elseif strcmp(Exper, 'JS'), CxtList = {'Forest','','City'};
            if thisRip.area==5, thisRip.area=0; end
        end
        cxt = CxtList{thisRip.context};
        if thisRip.correctness, corr='Correct'; else, corr='Wrong'; end
        title(['trial ' (thisRip.trial{1}(end-2:end)) ', ' cxt ', ' corr],'fontsize',15)
        axis off
        
        % position
        subplot(9,6,6)
        title([jjnum2str(thisRip.speed,2) 'cm/s'],'fontsize',15)
        axis off
        %%
        % reactivated cells
        %         [~,m] = max(abs(RP1));
        subplot(9,6,[7,8,13,14])
        thisEEG = EEG.(['TT' num2str(TargetTT_p)]).Raw(Ist:Ied);
        plot(thisEEG,'k')
        %         hold on
        %         plot(EEG.(['TT' num2str(TargetTT_p)]).Filtered(Ist:Ied),'b')
        
        title(['TT' num2str(TargetTT_p)])
        x1=mar; x2=mar+thisRip.RippleDuration*Params.Fs;
        xlim([0 dur])
        line([x1 x1], [min(thisEEG) max(thisEEG)+50], 'color','k','linestyle','--')
        line([x2 x2], [min(thisEEG) max(thisEEG)+50], 'color','k','linestyle','--')
        axis off
                
        
        %%
        spks_epoch=[];u=0; Units={};UnitsA={};
        Clist=jet(256);
        cls = size(clusters,1);
        for un = 1:cls
            thisTTID = num2str(clusters.TT(un));
            thisCLID = num2str(str2double(clusters.ID{un}(end-1:end)));
            Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);
            
            thisSpks = Spk.t_spk(Spk.t_spk>=thisRip.STtime-mar/Params.Fs & Spk.t_spk<=thisRip.EDtime+mar/Params.Fs);
            thisSpks_in = Spk.t_spk(Spk.t_spk>=thisRip.STtime & Spk.t_spk<=thisRip.EDtime);
            if ~isempty(thisSpks_in)
                spks_epoch = [spks_epoch;[thisSpks,ones(size(thisSpks,1),1)*un,ones(size(thisSpks,1),1)*u,ones(size(thisSpks,1),1)*clusters.SI(un),...
                    ones(size(thisSpks,1),1)*clusters.RDI_LScene(un), ones(size(thisSpks,1),1)*clusters.RDI_RScene(un), ones(size(thisSpks,1),1)*clusters.RDI_LR(un)]];
                u=u+1;
                Units = [Units; [thisTTID '-' thisCLID]];
                UnitsA = [UnitsA; [clusters.ID(un)]];
            end
            
        end
        %%
        spks_epochA=[];u=0; 
        cls = size(clusters_A,1);
        
        for un = 1:cls
            thisTTID = num2str(clusters_A.TT(un));
            thisCLID = num2str(str2double(clusters_A.ID{un}(end-1:end)));
            Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);
            
            thisSpks = Spk.t_spk(Spk.t_spk>=thisRip.STtime-mar/Params.Fs & Spk.t_spk<=thisRip.EDtime+mar/Params.Fs);
            thisSpks_in = Spk.t_spk(Spk.t_spk>=thisRip.STtime & Spk.t_spk<=thisRip.EDtime);
            if ~isempty(thisSpks_in)

                spks_epochA = [spks_epochA;[thisSpks,ones(size(thisSpks,1),1)*un,ones(size(thisSpks,1),1)*u ,ones(size(thisSpks,1),1)*clusters_A.SI(un),...
                    ones(size(thisSpks,1),1)*clusters_A.RDI_LScene(un), ones(size(thisSpks,1),1)*clusters_A.RDI_RScene(un), ones(size(thisSpks,1),1)*clusters_A.RDI_LR(un)]];
                u=u+1;
            end
        end
        %%
        subplot(9,6,[19,20,25,26])
        hold on
        for s=1:size(spks_epoch,1)
            if ismember(spks_epoch(s,1),spks_epochA)
                if spks_epoch(s,5)>0, cl='r';
                elseif spks_epoch(s,5)<0, cl='b';
                    else, cl='k';
                end
            else, cl='k'; end

            x = (spks_epoch(s,1) - thisRip.STtime)*1e3+100;
            patch([x-ti x+ti x+ti x-ti], [spks_epoch(s,3)+.2 spks_epoch(s,3)+.2 spks_epoch(s,3)+.8 spks_epoch(s,3)+.8],...
                cl,'edgecolor',cl)
        end
        x1=mar/2; x2=thisRip.RippleDuration*1e3+mar/2;
        xlim([0 dur/2])
        line([x1 x1], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
        line([x2 x2], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
        text(x1+1,0,'0')
        text(x2+1,0,[num2str(thisRip.RippleDuration*1e3) 'ms'])
        axis off
         %%
        subplot(9,6,[31 32 37 38])
        hold on
        for s=1:size(spks_epoch,1)
            if ismember(spks_epoch(s,1),spks_epochA)
                if spks_epoch(s,6)>0, cl='r';
                elseif spks_epoch(s,6)<0, cl='b';
                else, cl='k';
                end
            else, cl='k'; end
            x = (spks_epoch(s,1) - thisRip.STtime)*1e3+100;
            patch([x-ti x+ti x+ti x-ti], [spks_epoch(s,3)+.2 spks_epoch(s,3)+.2 spks_epoch(s,3)+.8 spks_epoch(s,3)+.8],...
                cl,'edgecolor',cl)
        end

        x1=mar/2; x2=thisRip.RippleDuration*1e3+mar/2;
        xlim([0 dur/2])
        line([x1 x1], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
        line([x2 x2], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
        text(x1+1,0,'0')
        text(x2+1,0,[num2str(thisRip.RippleDuration*1e3) 'ms'])
        axis off
         %%
        subplot(9,6,[43 44 49 50])
        hold on
        for s=1:size(spks_epoch,1)
            if ismember(spks_epoch(s,1),spks_epochA)
                if spks_epoch(s,7)>0, cl='r';
                elseif spks_epoch(s,7)<0, cl='b';
                    else, cl='k';
                end
            else, cl='k'; end
            x = (spks_epoch(s,1) - thisRip.STtime)*1e3+100;
            patch([x-ti x+ti x+ti x-ti], [spks_epoch(s,3)+.2 spks_epoch(s,3)+.2 spks_epoch(s,3)+.8 spks_epoch(s,3)+.8],...
                cl,'edgecolor',cl)
        end

        x1=mar/2; x2=thisRip.RippleDuration*1e3+mar/2;
        xlim([0 dur/2])
        line([x1 x1], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
        line([x2 x2], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
        text(x1+1,0,'0')
        text(x2+1,0,[num2str(thisRip.RippleDuration*1e3) 'ms'])
        axis off
        %%
        subplot(9,6,[9 10 15 16])
        hold on
        for s=1:size(spks_epoch,1)
            x = (spks_epoch(s,1) - thisRip.STtime)*1e3+100;
            patch([x-ti x+ti x+ti x-ti], [spks_epoch(s,3)+.2 spks_epoch(s,3)+.2 spks_epoch(s,3)+.8 spks_epoch(s,3)+.8],...
                'k','edgecolor','k')
        end
        x1=mar/2; x2=thisRip.RippleDuration*1e3+mar/2;
        xlim([0 dur/2])
        line([x1 x1], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
        line([x2 x2], [min(spks_epoch(:,3)) max(spks_epoch(:,3))+1], 'color','k','linestyle','--')
        axis off
       [~,ia,~] = unique(spks_epoch(:,2:3),'rows');
       spks_epoch_u = spks_epoch(ia,:);
%        spks_epoch_u(isnan(spks_epoch_u))='';
        for cl=1:size(Units,1)
            text(400,-0.5+cl,[Units{cl} '    ' jjnum2str(spks_epoch_u(cl,5),2) '    ' jjnum2str(spks_epoch_u(cl,6),2) '    ' jjnum2str(spks_epoch_u(cl,7),2)])
        end
        
        %%
        
        FRMap = [];
        for cl=1:size(UnitsA)
            try
            thisField = table;
            temp=[];
            thisTTID = num2str(str2double(UnitsA{cl}(8:9)));
            thisCLID = num2str(str2double(UnitsA{cl}(11:12)));
            Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);
            temp = Spk.t_spk(~Spk.area_spk(:,5) & Spk.correctness_spk(:,1));
            temp = sortrows(temp,1);
            
            thisField.ts = temp(:,1);
            thisMap.thisFieldMap1D = getFieldMaps(UnitsA{cl},thisField,'session',ROOT.Raw.Mother,ROOT.Info);
            
            
            for i=1:numel(thisMap.thisFieldMap1D.skaggsMap1D)
                FRMap(i,:,cl) = thisMap.thisFieldMap1D.skaggsMap1D{i};
            end
            
            FRMap(6,:,cl) = thisMap.thisFieldMap1D.skaggsMap_left1D;
            FRMap(7,:,cl) = thisMap.thisFieldMap1D.skaggsMap_right1D;
            
            FRMap(:,:,cl) =1- FRMap(:,:,cl) / max(max(FRMap(:,:,cl)));
            catch
                FRMap(:,:,cl) =nan;
            end
        end
        FRMap(isnan(FRMap))=1;
        %%
        subplot(9,6,[21 22 27 28])
        imagesc(flip(squeeze(FRMap(2,:,:))'))
        
        caxis([0 max(max(max(FRMap)))])
        axis off
        for i=1:size(FRMap,3)
            text(-8,i,Units{size(FRMap,3)-i+1})
        end
        title(CxtList{1},'color','r')
        
        
        subplot(9,6,[23 24 29 30])
        imagesc(flip(squeeze(FRMap(3,:,:))'))
        caxis([0 max(max(max(FRMap)))])
        axis off
        axis off
        for i=1:size(FRMap,3)
            text(size(FRMap,2),i,jjnum2str(spks_epoch_u(size(FRMap,3)-i+1,5),2))
        end
        title(CxtList{3},'color','b')
        
        if ~strcmp(Exper,'JS')
            
            subplot(9,6,[33 34 39 40])
            imagesc(flip(squeeze(FRMap(4,:,:))'))
            caxis([0 max(max(max(FRMap)))])
            axis off
            for i=1:size(FRMap,3)
                text(-8,i,Units{size(FRMap,3)-i+1})
            end
            title('Pebbles','color','r')
            
            subplot(9,6,[35 36 41 42])
            imagesc(flip(squeeze(FRMap(5,:,:))'))
            caxis([0 max(max(max(FRMap)))])
            axis off
            for i=1:size(FRMap,3)
                text(size(FRMap,2),i,jjnum2str(spks_epoch_u(size(FRMap,3)-i+1,6),2))
            end
            title('Mountains','color','b')
            
            subplot(9,6,[45 46 51 52])
            imagesc(flip(squeeze(FRMap(6,:,:))'))
            caxis([0 max(max(max(FRMap)))])
            axis off
            for i=1:size(FRMap,3)
                text(-8,i,Units{size(FRMap,3)-i+1})
            end
            title('Left','color','r')
            
            subplot(9,6,[47 48 53 54])
            imagesc(flip(squeeze(FRMap(7,:,:))'))
            caxis([0 max(max(max(FRMap)))])
            axis off
            for i=1:size(FRMap,3)
               text(size(FRMap,2),i,jjnum2str(spks_epoch_u(size(FRMap,3)-i+1,7),2))
            end
            title('Right','color','b')
        end
        
        colormap(gray)
        %%

%         fd = [ROOT.Fig3 '\' ['rat' jmnum2str(thisRip.rat,3)] '\rat' thisRSID];
%         if ~exist(fd), mkdir(fd); end
        saveas(gca,[ROOT.Fig3 '\' cell2mat(thisRip.ID) '.png'])

        close all
    catch
        saveas(gca,[ROOT.Fig3 '\' cell2mat(thisRip.ID) '.png'])
        close all
    end
    
end
