function [thisRID,thisSID,thisTTID,thisCLID, thisFLID] = parsing_clusterID(clusterID,d)

hyp = find(clusterID=='-');

thisRID = clusterID(1:hyp(1)-1);
thisSID=0;
thisTTID=0;
thisCLID = 0;
thisFLID = 0;
switch length(hyp)
    case 1
thisSID = clusterID(hyp(1)+1:end);
    case 2
    thisSID= clusterID(hyp(1)+1:hyp(2)-1);
    thisTTID =  clusterID(hyp(2)+1:end);
    case 3
    thisSID= clusterID(hyp(1)+1:hyp(2)-1);
    thisTTID =  clusterID(hyp(2)+1:hyp(3)-1);
    thisCLID =  clusterID(hyp(3)+1:end);
    case 4
    thisSID= clusterID(hyp(1)+1:hyp(2)-1);
    thisTTID =  clusterID(hyp(2)+1:hyp(3)-1);
    thisCLID =  clusterID(hyp(3)+1:hyp(4)-1);
    thisFLID = clusterID(hyp(4)+1:end);
    otherwise
end

if d~=0
thisRID  = str2double(thisRID);
thisSID  = str2double(thisSID);
thisTTID  = str2double(thisTTID);
thisCLID  = str2double(thisCLID);
thisFLID  = str2double(thisCLID);

if d==1
thisRID  = num2str(thisRID);
thisSID  = num2str(thisSID);
thisTTID  = num2str(thisTTID);
thisCLID  = num2str(thisCLID);
thisFLID  = num2str(thisCLID);
elseif d==2
thisRID  =jmnum2str(thisRID,3);
thisSID  = jmnum2str(thisSID,2);
thisTTID  = jmnum2str(thisTTID,2);
thisCLID  = jmnum2str(thisCLID,2);
thisFLID  = jmnum2str(thisCLID,2);
end
end
