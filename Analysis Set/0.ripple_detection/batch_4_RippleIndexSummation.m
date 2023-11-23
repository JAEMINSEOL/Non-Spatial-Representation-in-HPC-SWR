
Initial_SWRFilter_common;
warning off
thisRegion = 'SUB';
Experimenter = {'JS','LSM','SEB'};

ROOT.Save = ROOT.Processed;
ROOT.Old = [ROOT.Save '\ripples_mat\R0'];
ROOT.Unit = [ROOT.Save '\units_mat\U0'];
ROOT.SaveRip = [ROOT.Save '\ripples_mat\R0'];
ROOT.Fig = [ROOT.Save '\ripples_mat\ProfilingSheet\R0\' thisRegion];
if ~exist(ROOT.SaveRip), mkdir(ROOT.SaveRip); end
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end

%%
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_' thisRegion '_filtered.xlsx'],'ReadVariableNames',true);


fd = dir(ROOT.Save);

RipplesTable_all = table;
% RipplesTable_all = readtable([ROOT.Info '\RipplesList_' thisRegion '.xlsx']);
%%
for sid=1:size(SessionList,1)
    if SessionList.include(sid) & ismember(SessionList.experimenter{sid},Experimenter)
            
        ID = [jmnum2str(SessionList.rat(sid),3) '-' jmnum2str(SessionList.session(sid),2)];
        Recording_region_TT = Recording_region({ID},:);
        thisRID = SessionList.rat(sid);
        thisSID = SessionList.session(sid);
        thisRSID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];
        
        TargetTT_r = GetTargetTT(ROOT,thisRSID,thisRegion,Params,0.03);
        TargetTT_c = GetTargetTT(ROOT,thisRSID,thisRegion,Params,0);
        if size(TargetTT_c,1)<2, continue; end
        if ~exist([ROOT.Old '\' thisRegion '\' ID '.mat']), continue; end
        load([ROOT.Old '\' thisRegion '\' ID '.mat']);
         %%
        L=1;
        for t1=1:length(TargetTT_c)
            try
                thisTTID=TargetTT_r(t1);
                L = max(L,length(Ripples.(['TT' num2str(thisTTID)]).index(end,2)));
                cscID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2) '-' num2str(thisTTID)];
                % load CSC data
                cscData = loadCSC(cscID, ROOT.Raw.Mother, Params.CSCfileTag, Params.exportMODE, Params.behExtraction,'CSC');
                Timestamps_expand = cscData.Timestamps ./ 10^6; % make to sec unit
                break;
            end
        end
   
        time_ori = Timestamps_expand(1);
%%         
        RippleArray = zeros(L,24);

        for t1=1:length(TargetTT_r)
            try
                t=TargetTT_r(t1);
                RippleArray(1:Ripples.(['TT' num2str(t)]).index(end,2),t) =...
                    jmIndex2Bin(Ripples.(['TT' num2str(t)]).index,Ripples.(['TT' num2str(t)]).index(end,2));
            end
        end
        RippleVector = sum(RippleArray,2);
%         %%
%         
%         n=1;ripples_index=[];
%         Params.NumTTthreshold=min(3,(round(length(TargetTT_r)/2)));
%         ripples_index = jmBin2Index(RippleVector,Params.NumTTthreshold);
%         
%         if isempty(ripples_index), disp([ID ' has no ripple!']); continue; end
%         
%         % duration check
% %         ripples_index(ripples_index(:,2)-ripples_index(:,1) < Params_Ripple.minDuration*Params_Ripple.Fs, :) = [];
%         
%         % grouping
%         ripples2=[]; ripples_index2=[];
%         i=1;
%         while i<size(ripples_index,1)
%             j=i+1;
%             while ripples_index(j,1) - ripples_index(j-1,2) <= Params_Ripple.groupingInterval*Params_Ripple.Fs
%                 j=j+1;
%                 if j>size(ripples_index,1), break; end
%             end
%             ripples_index2 = [ripples_index2;[ripples_index(i,1),ripples_index(j-1,2)]];
%             i=j;
%         end
%         if ripples_index(end,2)~=ripples_index2(end,2)
%             ripples_index2 = [ripples_index2 ; ripples_index(end,1:2)];
%         end
%         ripples_index= ripples_index2;
%         
        %%
        TT_table = readtable([ROOT.Info '\TT_table.xlsx']);
         thisTT_table = TT_table(TT_table.rat==thisRID & TT_table.session==thisSID& strcmp(TT_table.region,thisRegion),:);
         for t=1:size(thisTT_table,1)
             if ~ismember(thisTT_table.TT(t),TargetTT_r)
                 thisTT_table.TT(t)=0;
             end
         end
         thisTT_table= thisTT_table(thisTT_table.TT~=0,:);
          [~,t] = max(thisTT_table.RippleBandMean);
         TargetTT_p = thisTT_table.TT(t);
        
         ripples_index = Ripples.(['TT' num2str(TargetTT_p)]).index;
         %%
         % duration check
         ripples_index(ripples_index(:,2)-ripples_index(:,1) < Params_Ripple.minDuration*Params_Ripple.Fs, :) = [];
                ripples_index(ripples_index(:,2)-ripples_index(:,1) > Params_Ripple.maxDuration*Params_Ripple.Fs, :) = [];
         
        if ripples_index(end,2)<ripples_index(end,1), ripples_index(end,2) = size(RippleVector,1); end
         ripples_index(ripples_index==0) = 1;
%%
        ripples_index(:,3) = time_ori + ripples_index(:,1)/Params_Ripple.Fs;
        ripples_index(:,4) = time_ori + ripples_index(:,2)/Params_Ripple.Fs;
        for i=1:size(ripples_index,1), ripples_index(i,5) = max(RippleVector(ripples_index(i,1):ripples_index(i,2))); end
        
        writematrix(ripples_index,[ROOT.SaveRip '\' ID '_' thisRegion '.csv'],'WriteMode', 'overwrite')
        disp([ID '_' thisRegion '.csv is saved!'])
        
        
        ripples = readtable([ROOT.SaveRip '\' ID '_' thisRegion '.csv']);
        Spike=LoadSpikeData(ROOT, ID, TargetTT_c, Params.cellfindn);
        clusters = Cluster_List(Cluster_List.rat==thisRID & Cluster_List.session==thisSID,:);
        %%
        if isempty(clusters)
            spks_all=0;
        else
        spks_all=[];
        for cl = 1:size(clusters,1)
            thisTTID = num2str(clusters.TT(cl));
            thisCLID = num2str(str2double(clusters.ID{cl}(end-1:end)));
            Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);
            spks_all = [spks_all;[Spk.t_spk,ones(size(Spk.t_spk,1),1)*cl]];
        end
        end
                Ripple_List=table;
        
        for rid = 1: size(ripples,1)
            Ripple_List.ID{rid} = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2) '-' thisRegion '-' jmnum2str(rid,4)];
            Ripple_List.rat(rid) = thisRID;
            Ripple_List.session(rid) = thisSID;
            Ripple_List.experimenter{rid} = SessionList.experimenter{sid};
            Ripple_List.region{rid} = thisRegion;
            Ripple_List.ripple(rid) = rid;
            
            Ripple_List.NumAllTT(rid) = size(TargetTT_r,1);
            temp_r = RippleArray(ripples.Var1(rid):ripples.Var2(rid),:);
           [m,t]  = max(RippleVector(ripples.Var1(rid):ripples.Var2(rid)));
             Ripple_List.NumTT(rid) = m;
            Ripple_List.TTs(rid) = binaryVectorToDecimal(flip(temp_r(t,:)));
      

            Ripple_List.STindex(rid) = ripples.Var1(rid);
            Ripple_List.EDindex(rid) = ripples.Var2(rid);
            Ripple_List.STtime(rid) = ripples.Var3(rid);
            Ripple_List.EDtime(rid) = ripples.Var4(rid);
            
            spks = find((spks_all(:,1)>=Ripple_List.STtime(rid)) & (spks_all(:,1)<=Ripple_List.EDtime(rid)));
            Ripple_List.spike(rid) = length(spks);
            if ~isempty(spks)
                Ripple_List.ensemble(rid) = length(unique(spks_all(spks,2)));
            else
                Ripple_List.ensemble(rid) = 0;
            end
            
        end
        
        RipplesTable_all = [RipplesTable_all;Ripple_List];
        
        %%
        ripples.ensemble = Ripple_List.ensemble;
        ID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];
%         EEG = LoadEEGData(ROOT, ID, TargetTT_r,Params,Params_Ripple);
%         RipplePropertySheet_Lv1_exc(ROOT,SessionList(sid,:),ripples,EEG,TargetTT_r,Spike,clusters,Params_Ripple);
   
        disp([ID ' profile sheet is finished!'])
    end
end
%%
writetable(RipplesTable_all,[ROOT.SaveRip '\RipplesList_' thisRegion '.xlsx'],'WriteMode', 'overwrite')
RipplesTable_filt = RipplesTable_all(RipplesTable_all.ensemble>=3,:);
writetable(RipplesTable_filt,[ROOT.SaveRip '\RipplesList_' thisRegion '_filtered.xlsx'],'WriteMode', 'overwrite')
