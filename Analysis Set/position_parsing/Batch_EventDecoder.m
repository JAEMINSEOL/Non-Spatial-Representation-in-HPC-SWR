
Initial_SWRFilter_common;
warning off

Session_List = readtable([ROOT.Info '\SessionList_SWR.xlsx']);

Epoch=zeros(15,2);
stCellRun = 1;

for sid = 1 : size(Session_List,1)
    
    if (strcmp(Session_List.experimenter{sid},'LSM') || strcmp(Session_List.experimenter{sid},'SEB') ) && Session_List.include(sid)
        
        try
            thisRID = jmnum2str(Session_List.rat(sid),3);
            thisSID = jmnum2str(Session_List.session(sid),2);
            
%             EventDecoderFunction_forJS(thisRID,thisSID,ROOT);
            Behav = EditParsedPos(ROOT,Session_List,sid);
 
            disp([thisRID '-' thisSID 'is processed!'])
        catch
            disp([thisRID '-' thisSID 'is failed!'])
        end
%         Epoch = MkBehavEpoch(ROOT,Session_List,sid,Epoch);
        
        
    end
end