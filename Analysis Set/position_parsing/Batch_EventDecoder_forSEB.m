
clc;
clear;
fclose all;

%
mother_root = 'F:\EPhysRawData\RawData';
addpath(genpath('D:\HPC-SWR project\Analysis program'))

sessionID = {'052-01'};
% sessionID = {'232-04', '232-05', '232-06', '232-07', '234-01', '234-02', '234-03', '234-04'};
% sessionID = {'295-01','295-02','295-03','295-04'};
% sessionID = {'232-13','234-09','295-10','415-18'};
% sessionID = {'232-08','234-05','295-05','415-14','561-06'};
% ,...
%     '313-03','313-04','425-03','425-04','454-03','454-04','471-03','471-04',...
%     '487-03','487-04','553-03','553-04','562-03','562-04'};
end2Flag = 0;   % if there's no end2 signal, then end2Flag = 0. 
                % if there's end2 signal(3rd contexts pair) follwed by end1 signal, then end2Flag = 1.
                % if there's end2 signal(blurred or trace) follwed by end1 siganl, then end2Flag = 2.
track_type = 2; % 1 for short track, 2 for long track.
%

for ssRUN = 1 : length(sessionID)
    
    findHYPHEN = find(sessionID{ssRUN} == '-');
    thisRID = sessionID{ssRUN}(1, 1:findHYPHEN(1) - 1);
    thisSID = sessionID{ssRUN}(1, findHYPHEN(1) + 1:end);
    
    session_root = [mother_root '\rat' thisRID '\rat' thisRID '-' thisSID];    
    EventDecoderFunction(session_root,end2Flag, track_type);
    CheckEventParsingFunction(session_root, thisRID, thisSID);
end

posFileName = 'VT1.nvt';