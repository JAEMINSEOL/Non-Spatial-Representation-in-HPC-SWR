% WriteBehaviorEpoch.csv using Parsed Position


Initial_SWRFilter_common;

% mother_root = 'H:\CA1&SUB_SCSM\ephys_analysis';
mother_root = ['F:\EPhysRawData\RawData'];

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv']);

Cluster_list = table;
exper = {'SEB'};
TargRegion = 'CA1';
Epoch=zeros(15,2);
for sid = 1:size(Session_List,1)
    temp = table;
    if ismember(Session_List.experimenter{sid},exper)
        clear t
        thisRID = Session_List.rat(sid);
        thisSID = Session_List.session(sid);
        try
            loc = [ROOT.Raw.Mother '\rat' jmnum2str(thisRID,3) '\rat' jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];
            load([loc '\Behavior\ParsedPosition.mat'])
            d = dir([loc '\Behavior\ParsedPosition.mat']);
            cd([loc '\Behavior'])
            copyfile(d.name, [loc '\ParsedPosition.mat'])
            Epoch(thisSID,1) = t(1)*10^6;
            Epoch(thisSID,2) = t(end)*10^6;
        catch
        end
        
        
        if Session_List.rat(sid+1)~=thisRID
            xlswrite([ROOT.Raw.Mother '\rat' jmnum2str(thisRID,3) '\behaviorEpoch_rat' jmnum2str(thisRID,3) '.xlsx'],Epoch)
            Epoch=zeros(15,2);
        end
    end
end
        