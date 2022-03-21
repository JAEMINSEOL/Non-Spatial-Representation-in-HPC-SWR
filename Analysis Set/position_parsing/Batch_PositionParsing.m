% Batch_PositionParsing

clc;
clear;
fclose all;

%
mother_root = 'D:\HPC-LFP project\RawData';

% sessionID = {'561-01','561-02','561-03','561-04','561-05','561-06'};
% sessionID = {'232-04', '232-05', '232-06', '232-07', '232-08', ...
%     '234-01', '234-02', '234-03', '234-04', '234-05'};
% sessionID = {'232-13','234-09','295-10','415-18'};
sessionID = {'561-01','561-02','561-03','561-04','561-05','561-06'};
%

for ssRUN = 1 : length(sessionID)
    
    findHYPHEN = find(sessionID{ssRUN} == '-');
    thisRID = sessionID{ssRUN}(1, 1:findHYPHEN(1) - 1);
    thisSID = sessionID{ssRUN}(1, findHYPHEN(1) + 1:end);
    
    cd([mother_root '\rat' thisRID]);
%     if str2double(thisRID) > 400
%         thisEPOCH = ['behaviorEpoch_rat' thisRID '.xlsx'];
%         epochST = xlsread(thisEPOCH, ['A' thisSID ':' 'A' thisSID]);
%         epochED = xlsread(thisEPOCH, ['B' thisSID ':' 'B' thisSID]);
%     else
        thisEPOCH = ['behaviorEpoch_rat' thisRID '.csv'];
        epochST = csvread(thisEPOCH, str2double(thisSID)-1, 0, [str2double(thisSID)-1, 0, str2double(thisSID)-1, 0]);
        epochED = csvread(thisEPOCH, str2double(thisSID)-1, 1, [str2double(thisSID)-1, 1, str2double(thisSID)-1, 1]);
%     end
%     
    session_root = [mother_root '\rat' thisRID '\rat' thisRID '-' thisSID];
    cd(session_root);
    
    createParsedPosition(session_root, epochST, epochED);
end