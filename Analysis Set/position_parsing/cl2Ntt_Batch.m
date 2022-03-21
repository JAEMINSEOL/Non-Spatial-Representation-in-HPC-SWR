
clear; clc; fclose all;

%
% mother_drive = 'D:\Hyunwoo Lee\backup(H)\CA1&SUB_SCSM';
mother_drive = 'D:\SUB-CA1 Ephys';

addpath([mother_drive '\1. Analysis program\SPK analysis programs\2. position spike parsing(final)']);

% mother_root = 'H:\CA1&SUB_SCSM\ephys_analysis';
mother_root = [mother_drive '\2. SPK Data'];

%
cd(mother_root);
[~,inputCSV] = xlsread('temp_rat415.csv');

% outputCSV = fopen('clusterID.csv', 'a');

stCellRun = 1;

for cellRUN = stCellRun : 1 : size(inputCSV,1)
        
    clusterID = inputCSV{cellRUN};

%     findHYPHEN = find(thisCLUSTER == '-');    
%     thisRID = (thisCLUSTER(1, 1:findHYPHEN(1) - 1));
%     thisSID = (thisCLUSTER(1, findHYPHEN(1) + 1:findHYPHEN(2) - 1));
%     thisTTID = (thisCLUSTER(1, findHYPHEN(2) + 1:findHYPHEN(3) - 1));
%     thisCLID = (thisCLUSTER(1, findHYPHEN(3) + 1:end));
%     
%     if thisSID < 10
%         clusterID = [num2str(thisRID) '-0' num2str(thisSID) '-' num2str(thisTTID)];
%     else
%         clusterID = [num2str(thisRID) '-' num2str(thisSID) '-' num2str(thisTTID)];
%     end
%     
%     if thisCLID < 10
%         clusterID = [clusterID '-0' num2str(thisCLID)];
%     else
%         clusterID = [clusterID '-' num2str(thisCLID)];
%     end
    
    cl2Ntt(mother_root, clusterID, clusterID);

%     fprintf(outputCSV, [clusterID '\n']);
    
    disp([clusterID ' has processed.']);
end


fclose all;