Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R0'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx']);
thisRegion = 'CA1';
Experimenter = ['LSM''JS','SEB'];
RipplesTable_all = table;

for sid=1:size(SessionList,1)
    if SessionList.include(sid) & ismember(SessionList.experimenter{sid},Experimenter) 
        thisRID = SessionList.rat(sid);
        thisSID = SessionList.session(sid);
        ID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];
        
        Recording_region_TT = Recording_region({ID},:);
        TargetTT = find(cellfun(cellfind(thisRegion),table2array(Recording_region_TT)'));
        
        ripples = readtable([ROOT.Old '\' ID '_' thisRegion '.csv']);
        EEG = LoadEEGData(ROOT, ID, TargetTT,Params,Params_Ripple);
            Spike=LoadSpikeData(ROOT, ID, TargetTT, cellfindn);
            
        clusters = Cluster_List(Cluster_List.rat==thisRID & Cluster_List.session==thisSID,:);
       sheet = RipplePropertySheet_Lv1_exc(ROOT,ID,ripples,EEG,TargetTT,Spike,clusters);
    end
end
        