Initial_SWRFilter_common;

% mother_root = 'H:\CA1&SUB_SCSM\ephys_analysis';
mother_root = ['F:\EPhysRawData\RawData'];

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);

Cluster_list = table;
exper = 'SEB';

for sid = 175:size(Session_List,1)
    temp = table;
    if Session_List.include(sid) && strcmp(Session_List.experimenter(sid),exper)
        temp.rat = Session_List.rat(sid);
        temp.session = Session_List.session(sid);
        for thisTTID = 1:24
            temp.TT = thisTTID;
            fd = dir([ROOT.Raw.Mother '\rat' jmnum2str(temp.rat,3) '\rat' jmnum2str(temp.rat,3) '-' jmnum2str(temp.session,2) '\' 'TT' num2str(thisTTID)]);
            for fid = 1: size(fd,1)
                name = fd(fid).name;
                if contains(name,'_cluster')
                    try
                    findDOT = find(name == '.');
                    temp.unit = str2double(name(findDOT+1:end));
                    temp.experimenter = exper;
                    temp.ID = [jmnum2str(temp.rat,3) '-' jmnum2str(temp.session,2) '-' jmnum2str(thisTTID,2) '-' jmnum2str(temp.unit,2)];
                    temp = movevars(temp, 'ID', 'Before', 'rat');
                              
                    clusterID = [jmnum2str(temp.rat,3) '-' jmnum2str(temp.session,2) '-' num2str(thisTTID) '-' num2str(temp.unit)];
                     cl2Ntt(mother_root, clusterID, clusterID);
                     
                     Cluster_list = [Cluster_list; temp];
                    catch
                    end
                end
            end
        end
    end
end
writetable(Cluster_list,[ROOT.Info '\ClusterList_SWR.xlsx'])      