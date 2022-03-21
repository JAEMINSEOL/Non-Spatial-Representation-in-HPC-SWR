function [position, bigTimeStamp] = simplePosPars(posFile, EpochST, EpochED)



%% set only behavior epochn using default file (u have to use winclust)
% this part is not used anymore

% root = 'D:\CA3&CA1&Sub_recording_CCSC\cheetah data for Parsing\Rat234';
% defaultFile = [root '\' 'behaviorEpoch_Rat234.csv']; 
% 
% sessionName = 'Rat234-01';
% sessionNumber = str2num(sessionName(8:9));
% 
% row = sessionNumber - 1;
% 
% epochST = csvread(defaultFile,row,0,[row,0,row,0]); % Epoch start
% epochED = csvread(defaultFile,row,1,[row,1,row,1]); % Epoch end


%% load position in behavior epoch

fieldSelectionFlags = [1 1 1 1 0 0];
% 1: timestamp, 2: x position, 3: y position, 4: head direction, 5: targets, 6: points

% [Timestamp,xPos,yPos,HeadDirection] = position_down_sampling(posFile, EpochST, EpochED);

[Timestamp,xPos,yPos,HeadDirection] = Nlx2MatVT(posFile, fieldSelectionFlags, 0, 4, [EpochST EpochED]);

[Timestamp,xPos,yPos,HeadDirection] = position_down_sampling(Timestamp, xPos, yPos, HeadDirection);

% 3rd input: 0 for no header extraction / 1 for header extraction
% 4th input: extraction range mode / 4 for timestamp range

position = struct( ...
    'timeStamp',0, ...
    'xPos', 0, ...
    'yPos', 0, ...
    'dir', 0);

RecordIndex = length(Timestamp);

for iter = 1:RecordIndex
    position(iter).timeStamp = Timestamp(1,iter);
    position(iter).xPos = xPos(1,iter);
    position(iter).yPos = yPos(1,iter);
    position(iter).dir = HeadDirection(1,iter);
end

bigTimeStamp = position(end).timeStamp;



%% SB's original program
% I only left some useful codes

% if(fd < 3)
%     disp('error : please check posision file');
% end


% position = struct( ...
%     'timeStamp',0, ...
%     'xPos', 0, ...
%     'yPos', 0, ...
%     'dir', 0);
% 
%     position(iter).timeStamp = first;
%     position(iter).xPos = fread(fd, 1, 'float');
%     position(iter).yPos = fread(fd, 1, 'float');
%     position(iter).dir = fread(fd, 1, 'float');
%   
% bigTimeStamp = position(iter - 1).timeStamp;