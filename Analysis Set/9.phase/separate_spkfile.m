%% separate clusters vis theta phase

warning off
Initial_SWRFilter_common;

ROOT.phase = 'X:\E-Phys Analysis\HPC-LFP project\backup\Theta phase\mat files (Theta phase)(2020.12.23)';
Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);
TargRegion = 'SUB';
exper = {'LSM'};

Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_' TargRegion '.xlsx']);

for cid = 1:size(Cluster_List,1)
    if ismember(Cluster_List.experimenter{cid}, exper)
        try
        clusterID = Cluster_List.ID{cid};
[thisRID,thisSID,thisTTID,thisCLID] = parsing_clusterID(clusterID,1);
thisSID = jmnum2str(str2double(thisSID),2);
load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' thisCLID '.mat'])
load([ROOT.phase '\rat' thisRID '-' thisSID '-' thisTTID '-' jmnum2str(str2double(thisCLID),2) '.mat'])


cluster_phase=zeros(size(t_spk,1),1);
phase_spk = nan(size(t_spk,1),1);
for t=1:size(t_spk,1)
    c = find(t_spk(t)==thisPHASE.ts);
    if ~isempty(c)
        cluster_phase(t,1) = thisPHASE.cluster(c(1));
    end

    c=find(t_spk(t)==PHASE_mat.ts);
    if ~isempty(c)
    phase_spk(t) = PHASE_mat.phase(c(1));
end
end

save([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' thisCLID '.mat'], 'phase_spk', 'cluster_phase','-append');

for pid = 1:max(cluster_phase)
id = find(cluster_phase==pid);
phase_field = struct;
phase_field.thisFieldMap = thisFieldMap{pid,1};
phase_field.x_spk = x_spk(id,:);
phase_field.y_spk = y_spk(id,:);
phase_field.t_spk = t_spk(id,:);
phase_field.a_spk = a_spk(id,:);
phase_field.area_spk = area_spk(id,:);
phase_field.ambiguity_spk = ambiguity_spk(id,:);
phase_field.correctness_spk = correctness_spk(id,:);
phase_field.cont_spk = cont_spk(id,:);
phase_field.trial_spk = trial_spk(id,:);
phase_field.side_spk = side_spk(id,:);
phase_field.phase_spk = phase_spk(id,:);

save([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' thisCLID '_f' num2str(pid) '.mat'],'-struct', 'phase_field');
end
        catch
            disp([clusterID ' is failed!'])
        end
    end
end