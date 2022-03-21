
function cl2Ntt(mother_root, winID, ss3dID)

% load data
findHYPHEN = find(winID == '-');

thisRID = winID(1, 1:findHYPHEN(1) - 1);
thisSID = winID(1, findHYPHEN(1) + 1:findHYPHEN(2) - 1);
thisTTID = winID(1, findHYPHEN(2) + 1:findHYPHEN(3) - 1);
thisCLID = winID(1, findHYPHEN(3) + 1:end);

winID_root = [mother_root '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID];
cd(winID_root);

thisEpochCLTS = dlmread(['TT' thisTTID '_cluster.' num2str(str2num(thisCLID))],',',13,0);
thisEpochCLTS = thisEpochCLTS(:,18);
thisEpochCLTS = transpose(thisEpochCLTS);

if exist(['TT',thisTTID,'_beh.ntt'], 'file')
    [wholeCLTS, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(['TT',thisTTID,'_beh.ntt'],[1 1 1 1 1], 1, 1, 0);
elseif exist(['TT',thisTTID,'.ntt'], 'file')
    [wholeCLTS, ScNumbers, CellNumbers, Features, Samples, Header] = Nlx2MatSpike(['TT',thisTTID,'.ntt'],[1 1 1 1 1], 1, 1, 0);
end
% 

% select spikes
clusterINDEX = [];
for clRUN = 1:size(thisEpochCLTS, 2)
    clusterINDEX(clRUN) = find(thisEpochCLTS(clRUN) == wholeCLTS);
end

ScNumbers = ScNumbers(clusterINDEX);
CellNumbers = CellNumbers(clusterINDEX);
Features = Features(:, clusterINDEX);
Samples = Samples(:, :, clusterINDEX);
% 

% make new ntt file

findHYPHEN = find(ss3dID == '-');

thisRID = ss3dID(1, 1:findHYPHEN(1) - 1);
thisSID = ss3dID(1, findHYPHEN(1) + 1:findHYPHEN(2) - 1);
thisTTID = ss3dID(1, findHYPHEN(2) + 1:findHYPHEN(3) - 1);
thisCLID = ss3dID(1, findHYPHEN(3) + 1:end);

ss3dID_root = [mother_root '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID];
cd(ss3dID_root);

CellNumbers(1:end) = str2double(thisCLID);

Mat2NlxSpike(['TT' thisTTID '_beh_SS_' thisCLID '.ntt'], 0, 1, 0, [1 1 1 1 1 1], ...
    thisEpochCLTS, ScNumbers, CellNumbers, Features, Samples, Header);
% 


end