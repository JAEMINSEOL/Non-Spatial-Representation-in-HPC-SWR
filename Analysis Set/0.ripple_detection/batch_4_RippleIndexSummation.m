
Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R0_5'];
%%
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx']);
thisRegion = 'CA1';
Experimenter = {'LSM','SEB','JS'};



fd = dir(ROOT.Save);

RipplesTable_all = table;
RipplesTable_all = readtable([ROOT.Info '\RipplesList_' thisRegion '.xlsx']);
%%
for sid=32:size(SessionList,1)
    if SessionList.include(sid) & ismember(SessionList.experimenter{sid},Experimenter) 
        ID = [jmnum2str(SessionList.rat(sid),3) '-' jmnum2str(SessionList.session(sid),2)];
    Recording_region_TT = Recording_region({ID},:);
        thisRID = SessionList.rat(sid);
    thisSID = SessionList.session(sid);
    
    TargetTT = find(cellfun(cellfindn2(thisRegion),table2array(Recording_region_TT)'));
    if size(TargetTT,1)<2, continue; end
    load([ROOT.Old '\' ID '.mat']);
    %%
    for t1=1:length(TargetTT)
        try
        thisTTID=TargetTT(t1);
          L = length(Ripples.(['TT' num2str(thisTTID)]).index);
            cscID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2) '-' num2str(thisTTID)];
            % load CSC data
            cscData = loadCSC(cscID, ROOT.Raw.Mother, Params.CSCfileTag, Params.exportMODE, Params.behExtraction,'CSC');
            Timestamps_expand = cscData.Timestamps ./ 10^6; % make to sec unit
            break;
        end
    end

    time_ori = Timestamps_expand(1);
    
    RippleArray = zeros(L,size(TargetTT,1));
    %%
    for t1=1:length(TargetTT)
        try
        t=TargetTT(t1);
            RippleArray(1:length(Ripples.(['TT' num2str(t)]).index),t1) = Ripples.(['TT' num2str(t)]).index;
        end
    end
    RippleVector = sum(RippleArray,2);
    %%
    
    n=1;ripples_index=[]; Params.NumTTthreshold=max(3,floor(length(TargetTT)/2));
    for i=2:size(RippleVector,1)
        if RippleVector(i-1)<Params.NumTTthreshold && RippleVector(i)>=Params.NumTTthreshold
            ripples_index(n,1)=i;
        elseif RippleVector(i-1)>=Params.NumTTthreshold && RippleVector(i)<Params.NumTTthreshold
            ripples_index(n,2)=i;
            n=n+1;
        end
    end
    
    if isempty(ripples_index), disp([ID ' has no ripple!']); continue; end
    
        % duration check
    ripples_index(ripples_index(:,2)-ripples_index(:,1) < Params_Ripple.minDuration*Params_Ripple.Fs, :) = [];
    
    % grouping
    ripples2=[]; ripples_index2=[];
    i=1;
    while i<size(ripples_index,1)
        j=i+1;
        while ripples_index(j,1) - ripples_index(j-1,2) <= Params_Ripple.groupingInterval*Params_Ripple.Fs
            j=j+1;
            if j>size(ripples_index,1), break; end
        end
        ripples_index2 = [ripples_index2;[ripples_index(i,1),ripples_index(j-1,2)]];
        i=j;
    end
    if ripples_index(end,2)~=ripples_index2(end,2)
        ripples_index2 = [ripples_index2 ; ripples_index(end,1:2)];
    end
    ripples_index= ripples_index2;
    
        % duration check
    ripples_index(ripples_index(:,2)-ripples_index(:,1) < Params_Ripple.minDuration*Params_Ripple.Fs*2, :) = [];
    if ripples_index(end,2)<ripples_index(end,1), ripples_index(end,2) = size(RippleVector,1); end
    ripples_index(ripples_index==0) = 1;
    
    %%
    ripples_index(:,3) = time_ori + ripples_index(:,1)/Params_Ripple.Fs;
    ripples_index(:,4) = time_ori + ripples_index(:,2)/Params_Ripple.Fs;
    for i=1:size(ripples_index,1), ripples_index(i,5) = max(RippleVector(ripples_index(i,1):ripples_index(i,2))); end
    
%     writematrix(ripples_index,[ROOT.Old '\' ID '_' thisRegion '.csv'],'WriteMode', 'overwrite')
%     disp([ID '_' thisRegion '.csv is saved!'])

    
    ripples = readtable([ROOT.Old '\' ID '_' thisRegion '.csv']);
    Spike=LoadSpikeData(ROOT, ID, TargetTT, cellfindn);
    clusters = Cluster_List(Cluster_List.rat==thisRID & Cluster_List.session==thisSID,:);
    %%
    Ripple_List=table;
    spks_all=[];
    for cl = 1:size(clusters,1)
        thisTTID = num2str(clusters.TT(cl));
        thisCLID = num2str(str2double(clusters.ID{cl}(end-1:end)));
        Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);
        spks_all = [spks_all;[Spk.t_spk,ones(size(Spk.t_spk,1),1)*cl]];
    end
    
    for rid = 1: size(ripples,1)
        Ripple_List.ID{rid} = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2) '-' thisRegion '-' jmnum2str(rid,4)];
        Ripple_List.rat(rid) = thisRID;
        Ripple_List.session(rid) = thisSID;
        Ripple_List.region{rid} = thisRegion;
        Ripple_List.ripple(rid) = rid;
        
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
    EEG = LoadEEGData(ROOT, ID, TargetTT,Params,Params_Ripple);
    RipplePropertySheet_Lv1_exc(ROOT,SessionList(sid,:),ripples,EEG,TargetTT,Spike,clusters,Params_Ripple);
    
       disp([ID ' profile sheet is finished!'])
    end
end
%%
% writetable(RipplesTable_all,[ROOT.Info '\RipplesList_' thisRegion '.xlsx'],'writemode','overwrite')
% RipplesTable_all = RipplesTable_all(RipplesTable_all.ensemble>=3,:);
% writetable(RipplesTable_all,[ROOT.Info '\RipplesList_' thisRegion '_filtered.xlsx'],'writemode','overwrite')
