function [trialRecord, trialSize, header] = trialPars(summaryFile)

VOID = 6;

trialFieldN = 13;

fd = fopen(summaryFile, 'r');

if(fd < 3)
    disp('error : fail to open Event.csv file');
end

header = fgetl(fd);

trialRecordWithVoid = textscan(fd, '%d%d%d%d%d%d%f%f%f%f%f%f%f', 'Delimiter', ',');

trialSize = length(trialRecordWithVoid{1});

fclose(fd);

 trialRecordWithVoid{1,VOID}(trialRecordWithVoid{1,12}-trialRecordWithVoid{1,7}>6)=2; %trial timeout을 void 기준으로 추가(재민)

nonVoidTrialIter = 1;

for trialIter = 1:trialSize
    
    %delete void trial
    if trialRecordWithVoid{1,VOID}(trialIter,1) == 2
        continue;
    end
    
    %if it is not void trial
    for tFieldIter = 1:trialFieldN
        trialRecord{1,tFieldIter}(nonVoidTrialIter,1) = trialRecordWithVoid{1,tFieldIter}(trialIter,1);
    end
    nonVoidTrialIter = nonVoidTrialIter + 1;
end