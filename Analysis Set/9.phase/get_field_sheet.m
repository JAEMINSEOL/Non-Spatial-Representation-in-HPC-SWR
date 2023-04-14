%% separate clusters via theta phase

warning off
Initial_SWRFilter_common;

ROOT.phase = 'X:\E-Phys Analysis\HPC-LFP project\backup\mat files (phase mat)';
Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);
ClusterList_SM= readtable([ROOT.Info '\result csv field_property_FR-TP(master sheet).xlsx']);
TargRegion = 'SUB';
exper = {'LSM'};

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_' TargRegion '.xlsx']);
Cluster_List_p = Cluster_List;
Cluster_List_n = Cluster_List;
Cluster_List_Field = table;

for cid = 1:size(Cluster_List,1)
    if ismember(Cluster_List.experimenter{cid}, exper)
        try
            clusterID = Cluster_List.ID{cid};
            [thisRID,thisSID,thisTTID,thisCLID,~] = parsing_clusterID(clusterID,1);
            thisSID = jmnum2str(str2double(thisSID),2);

             id = find(strcmp([thisRID '-' thisSID '-' thisTTID '-' jmnum2str(str2double(thisCLID),2)],ClusterList_SM.clusterID)); 
                        if ~isempty(id),id=id(1); end

            fid = dir([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID]);
            numf=0;
            for f=1:size(fid,1)
                thisid = fid(f).name;
                if contains(thisid,['parsedSpike_' thisCLID '_'])
                    if ~isempty(find(thisid=='f'))
                        numf = str2double(thisid(find(thisid=='f')+1));
                    end
                end
            end
            Cluster_List_thisS=table;

            if numf>0
                for f = 1:numf
                    if f>numf, break; end
                    load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' thisCLID '_f' num2str(f) '.mat'])
                    Cluster_List_thisF =Cluster_List(cid,1:19);
                    Cluster_List_thisF.ID = [clusterID '-' num2str(f)];
%                     Cluster_List_thisF(1,2:end) = Cluster_List(cid,2:end);
                    Cluster_List_thisF.ReMap_LScene =nan;
                    Cluster_List_thisF.ReMap_RScene =nan;
                    Cluster_List_thisF.ReMap_LR =nan;

                    if isfield(thisFieldMap,'d')
                        if thisFieldMap.onmazeAvgFR1D(1)>thisFieldMap.onmazeAvgFR1D(2),...
                                Cluster_List_thisF.RDI_LScene = thisFieldMap.d(1); else, Cluster_List_thisF.RDI_LScene = -thisFieldMap.d(1); end

                        if thisFieldMap.onmazeAvgFR1D(3)>thisFieldMap.onmazeAvgFR1D(4),...
                                Cluster_List_thisF.RDI_RScene = thisFieldMap.d(2); else, Cluster_List_thisF.RDI_RScene = -thisFieldMap.d(2); end

                        if (thisFieldMap.onmazeAvgFR1D(1)+thisFieldMap.onmazeAvgFR1D(2))>(thisFieldMap.onmazeAvgFR1D(3)+thisFieldMap.onmazeAvgFR1D(4)),...
                                Cluster_List_thisF.RDI_LR = thisFieldMap.d(3); else, Cluster_List_thisF.RDI_LR = -thisFieldMap.d(3); end



                    else
                        if ~isempty(id)
                          if thisFieldMap.onmazeAvgFR1D(1)>thisFieldMap.onmazeAvgFR1D(2),...
                                Cluster_List_thisF.RDI_LScene = thisFieldMap.RMI(1); else, Cluster_List_thisF.RDI_LScene = -thisFieldMap.RMI(1); end

                        if thisFieldMap.onmazeAvgFR1D(3)>thisFieldMap.onmazeAvgFR1D(4),...
                                Cluster_List_thisF.RDI_RScene = thisFieldMap.RMI(2); else, Cluster_List_thisF.RDI_RScene = -thisFieldMap.RMI(2); end

                        if (thisFieldMap.onmazeAvgFR1D(1)+thisFieldMap.onmazeAvgFR1D(2))>(thisFieldMap.onmazeAvgFR1D(3)+thisFieldMap.onmazeAvgFR1D(4)),...
                                Cluster_List_thisF.RDI_LR = thisFieldMap.RMI(3); else, Cluster_List_thisF.RDI_LR = -thisFieldMap.RMI(3); end

                        else
                            Cluster_List_thisF.RDI_LScene=nan;
                            Cluster_List_thisF.RDI_RScene = nan;
                            Cluster_List_thisF.RDI_LR = nan;
                        end

                    end

                    if isfield(thisFieldMap,'COM') & numf>1

                        Cluster_List_thisF.PeakBin = thisFieldMap.COM(1);
                        Cluster_List_thisF.onMazeMaxFR_field = thisFieldMap.onmazeMaxFR1D(1);
                    else
                        if ~isempty(id)

                        Cluster_List_thisF.PeakBin = ClusterList_SM.COM(id);
                        Cluster_List_thisF.onMazeMaxFR_field = ClusterList_SM.peakFR(id);
                        Cluster_List_thisF.SI = ClusterList_SM.SIScore(id);

                        else
                            Cluster_List_thisF.PeakBin = nan;
                        Cluster_List_thisF.onMazeMaxFR_field = 0;
                        end


                    end

                    Cluster_List_Field = [Cluster_List_Field;Cluster_List_thisF];
                    Cluster_List_thisS = [Cluster_List_thisS; Cluster_List_thisF];
                    save([ROOT.Raw.Map '/rat' thisRID '-' thisSID '-' num2str(thisTTID) '-' jmnum2str(str2double(thisCLID),2) '-' jmnum2str(f,2) '.mat'],'-struct','thisFieldMap')
                    clear thisFieldMap
                end
            else
                try
                    load([ROOT.phase '\rat' thisRID '-' thisSID '-' thisTTID '-' jmnum2str(str2double(thisCLID),2) '.mat'])
                    thisFieldMap = thisFieldMap{1};
                    Cluster_List_thisF =Cluster_List(cid,1:19);
                    Cluster_List_thisF.ID = [clusterID '-1'];
                    Cluster_List_thisF.ReMap_LScene =nan;
                    Cluster_List_thisF.ReMap_RScene =nan;
                    Cluster_List_thisF.ReMap_LR =nan;
                    if thisFieldMap.onmazeAvgFR1D(2)>thisFieldMap.onmazeAvgFR1D(3),...
                            Cluster_List_thisF.RDI_LScene = thisFieldMap.d(1); else, Cluster_List_thisF.RDI_LScene = -thisFieldMap.d(1); end

                    if thisFieldMap.onmazeAvgFR1D(4)>thisFieldMap.onmazeAvgFR1D(5),...
                            Cluster_List_thisF.RDI_RScene = thisFieldMap.d(2); else, Cluster_List_thisF.RDI_RScene = -thisFieldMap.d(2); end

                    if (thisFieldMap.onmazeAvgFR1D(2)+thisFieldMap.onmazeAvgFR1D(3))>(thisFieldMap.onmazeAvgFR1D(4)+thisFieldMap.onmazeAvgFR1D(5)),...
                            Cluster_List_thisF.RDI_LR = thisFieldMap.d(3); else, Cluster_List_thisF.RDI_LR = -thisFieldMap.d(3); end



                    Cluster_List_thisF.PeakBin = thisFieldMap.COM(1);
                    Cluster_List_thisF.onMazeMaxFR_field = thisFieldMap.onmazeMaxFR1D(1);

                    save([ROOT.Raw.Map '/rat' thisRID '-' thisSID '-' num2str(thisTTID) '-' jmnum2str(str2double(thisCLID),2) '-' jmnum2str(1,2) '.mat'],'-struct','thisFieldMap')
                catch
                    %                     find(strcmp([thisRID '-' thisSID '-' num2str(thisTTID) '-' jmnum2str(str2double(thisCLID),2)], ClusterList_SM.clusterID))
                    Cluster_List_thisF =Cluster_List(cid,1:19);
                    Cluster_List_thisF.ID = [clusterID '-1'];
                    Cluster_List_thisF.ReMap_LScene =nan;
                    Cluster_List_thisF.ReMap_RScene =nan;
                    Cluster_List_thisF.ReMap_LR =nan;
                    Cluster_List_thisF.RDI_LScene =nan;
                    Cluster_List_thisF.RDI_RScene =nan;
                    Cluster_List_thisF.RDI_LR =nan;
                    Cluster_List_thisF.PeakBin = nan;
                    Cluster_List_thisF.onMazeMaxFR_field = 0;

                end
                Cluster_List_Field = [Cluster_List_Field;Cluster_List_thisF];
                Cluster_List_thisS = [Cluster_List_thisS; Cluster_List_thisF];
            end
            disp([clusterID ' is finished'])

        end
    end
%%
    
    Cluster_List.ReMap_LScene(cid) =nan;
    Cluster_List.ReMap_RScene(cid) =nan;
    Cluster_List.ReMap_LR(cid) =nan;

    [~,t] = max(abs(Cluster_List_thisS.RDI_LScene));
    Cluster_List.RDI_LScene(cid) = Cluster_List_thisS.RDI_LScene(t);

    [~,t] = max(abs(Cluster_List_thisS.RDI_RScene));
    Cluster_List.RDI_RScene(cid) = Cluster_List_thisS.RDI_RScene(t);

    [~,t] = max(abs(Cluster_List_thisS.RDI_LR));
    Cluster_List.RDI_LR(cid) = Cluster_List_thisS.RDI_LR(t);

    [~,t] = max(abs(Cluster_List_thisS.onMazeMaxFR_field));
    Cluster_List.PeakBin(cid) = Cluster_List_thisS.PeakBin(t);

    Cluster_List.NumField(cid) = size(Cluster_List_thisS,1);
%%
        Cluster_List_p.ReMap_LScene(cid) =nan;
    Cluster_List_p.ReMap_RScene(cid) =nan;
    Cluster_List_p.ReMap_LR(cid) =nan;

    [~,t] = max((Cluster_List_thisS.RDI_LScene));
    Cluster_List_p.RDI_LScene(cid) = Cluster_List_thisS.RDI_LScene(t);

    [~,t] = max((Cluster_List_thisS.RDI_RScene));
    Cluster_List_p.RDI_RScene(cid) = Cluster_List_thisS.RDI_RScene(t);

    [~,t] = max((Cluster_List_thisS.RDI_LR));
    Cluster_List_p.RDI_LR(cid) = Cluster_List_thisS.RDI_LR(t);

    [~,t] = max(abs(Cluster_List_thisS.onMazeMaxFR_field));
    Cluster_List_p.PeakBin(cid) = Cluster_List_thisS.PeakBin(t);

    Cluster_List_p.NumField(cid) = size(Cluster_List_thisS,1);
%%
        Cluster_List_n.ReMap_LScene(cid) =nan;
    Cluster_List_n.ReMap_RScene(cid) =nan;
    Cluster_List_n.ReMap_LR(cid) =nan;

    [~,t] = min((Cluster_List_thisS.RDI_LScene));
    Cluster_List_n.RDI_LScene(cid) = Cluster_List_thisS.RDI_LScene(t);

    [~,t] = min((Cluster_List_thisS.RDI_RScene));
    Cluster_List_n.RDI_RScene(cid) = Cluster_List_thisS.RDI_RScene(t);

    [~,t] = min((Cluster_List_thisS.RDI_LR));
    Cluster_List_n.RDI_LR(cid) = Cluster_List_thisS.RDI_LR(t);

    [~,t] = max(abs(Cluster_List_thisS.onMazeMaxFR_field));
    Cluster_List_n.PeakBin(cid) = Cluster_List_thisS.PeakBin(t);

    Cluster_List_n.NumField(cid) = size(Cluster_List_thisS,1);

end
Cluster_List(isnan(Cluster_List.RDI_LScene),:)=[];
Cluster_List_p(isnan(Cluster_List_p.RDI_LScene),:)=[];
Cluster_List_n(isnan(Cluster_List_n.RDI_LScene),:)=[];
Cluster_List_Field(isnan(Cluster_List_Field.RDI_LScene),:)=[];

writetable(Cluster_List,[ROOT.Save '\UnitsTable_' TargRegion '_forAnalysis.xlsx'],'writemode','overwrite');
writetable(Cluster_List_p,[ROOT.Save '\UnitsTable_' TargRegion '_p_forAnalysis.xlsx'],'writemode','overwrite');
writetable(Cluster_List_n,[ROOT.Save '\UnitsTable_' TargRegion '_n_forAnalysis.xlsx'],'writemode','overwrite');
writetable(Cluster_List_Field,[ROOT.Save '\UnitsTable_' TargRegion '_field_forAnalysis.xlsx'],'writemode','overwrite');
