% WriteBehaviorEpoch.csv using Parsed Position

function MkBehavEpoch(ROOT,Session_List,sid)
        clear t
        thisRID = Session_List.rat(sid);
        thisSID = Session_List.session(sid);
        try
            loc = [ROOT.Raw.Mother '\rat' jmnum2str(thisRID,3) '\rat' jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];
            load([loc '\ParsedPosition.mat'])
%             d = dir([loc '\Behavior\ParsedPosition.mat']);
%             cd([loc '\Behavior'])
%             copyfile(d.name, [loc '\ParsedPosition.mat'])
            Epoch(thisSID,1) = t(1);
            Epoch(thisSID,2) = t(end);
        catch
        end
        
        
        if Session_List.rat(sid+1)~=thisRID
            xlswrite([ROOT.Raw.Mother '\rat' jmnum2str(thisRID,3) '\behaviorEpoch_rat' jmnum2str(thisRID,3) '.xlsx'],Epoch)
            Epoch=zeros(15,2);
            disp(['behaviorEpoch_rat' jmnum2str(thisRID,3) '.xlsx is saved!'])
        end
    end
