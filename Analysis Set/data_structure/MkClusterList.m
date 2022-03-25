Initial_SWRFilter_common;

% mother_root = 'H:\CA1&SUB_SCSM\ephys_analysis';
mother_root = ['F:\EPhysRawData\RawData'];

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);

Cluster_list = table;
exper = {'SEB','LSM','JS'};
TargRegion = 'CA1';



for sid = 1:size(Session_List,1)
    temp = table;
    if Session_List.include(sid) && ismember(Session_List.experimenter{sid},exper)
        temp.rat = Session_List.rat(sid);
        temp.session = Session_List.session(sid);
        r = find(strcmp(Recording_region.SessionID,[jmnum2str(temp.rat,3) '-' jmnum2str(temp.session,2)]));
        for thisTTID = 1:24
            temp.TT = thisTTID;
            loc = [ROOT.Raw.Mother '\rat' jmnum2str(temp.rat,3) '\rat' jmnum2str(temp.rat,3) '-' jmnum2str(temp.session,2) '\' 'TT' num2str(thisTTID)];
            fd = dir(loc);
%             trans_model_forJS(loc,fd,thisTTID)
            fd = dir(loc);
            for fid = 1: size(fd,1)
                name = fd(fid).name;
                if contains(name,'_cluster') 
                    try
                    findDOT = find(name == '.');
                    temp.unit = str2double(name(findDOT+1:end));
                    temp.region{1} = Recording_region.(['TT' num2str(temp.TT)]){r};
                    temp.experimenter{1} = Session_List.experimenter{sid};
                    temp.date = Session_List.date(sid);
                    temp.session_type{1} = Session_List.type{sid};
                    temp.ID = [jmnum2str(temp.rat,3) '-' jmnum2str(temp.session,2) '-' jmnum2str(thisTTID,2) '-' jmnum2str(temp.unit,2)];
                    temp = movevars(temp, 'ID', 'Before', 'rat');
                    
                    clusterID = [jmnum2str(temp.rat,3) '-' jmnum2str(temp.session,2) '-' num2str(thisTTID) '-' num2str(temp.unit)];
                    %                       cl2Ntt(mother_root, clusterID, clusterID);

                    Cluster_list = [Cluster_list; temp];
                    catch
                    end
                end
            end
        end
    end
end
writetable(Cluster_list,[ROOT.Info '\ClusterList_SWR.xlsx'])
Cluster_list_CA1 = Cluster_list(strcmp(TargRegion,Cluster_list.region),:);
writetable(Cluster_list_CA1,[ROOT.Info '\ClusterList_SWR_' TargRegion '.xlsx'])

function trans_model_forJS(loc,fd,thisTTID)
            for fid = 1: size(fd,1)
                name = fd(fid).name;
                if (contains(name, ['T' num2str(thisTTID) '.']) && ~contains(name, 'TT')) && str2double(name(end))<100
                    cd(loc)
                d = dir([loc '\' name]);
                movefile(d.name, ['TT' num2str(thisTTID) '_cluster' name(end-1:end)])
                end
            end
end
