%%---------     main purpose of the program     -------------%%
%this function is to make ParsedPosition.mat file. this file contain
%variables below.

%revision 9/28/2013 by SB
%This program contains inbound data.
%'area' variable contains inbound position data in 4th column.
%ParsedPosition folder was not available in this revision. instead it use
%'ParsedPosition2' folder.

%     ------- example of result variables. -----------
% PositionCountMap [48*72*82] : double?. x??*y??*trial. ???? Pos2Map??? ???? ??.
% a [13669 * 1]: ? timestamp?? head direction ?? ??.
% correctness [13669 * 2]: ? position?? correctness ?? ?? [0 1] correct, [1 0] incorrect
% cont [13669 * 2]: ? position?? obj?? ?? obj1 : [1 0], obj2 : [0 1].
% originally it was obj. but I change it cont.
% side [12894 * 2]: position?? touch side ?? ?? 
% t [13669 * 1]: time stamp of each position
% total_trial_number [1*1]: total trial number
% trial [13669 * 82]: this matrix means which position is involved in each trials (position(ordered by time stamp) * trial number)
% trial_corr [82 * 1]: correctness of each trial
% trial_context [82 * 1]: object of each trial. original version, it was trial_context.
%But in this program object is context. So I changed it
% trial_void  [82 * 1]: void of trial. It was originally reversed.(1 * 82) I think this matrix reversed because
% of mistake. So I corrected it.
% x [13669 * 1]: x position of each position
% y [13669 * 1]: y position of each position

%                  sb added.
%trial_side [82 * 1] : touch side of each trial
%area[13669 * 3] %there are three area. It will be used with trial
%variables.

%made by SB
%%---------------    function prototype(if it is function)             --%%
function createParsedPosition(session_root, epochST, epochED)

% trial matrix notation.
INDEX = 1;
CORRECTNESS = 2;
CONTEXT = 3;
AMBIGUITY = 4;
RESPONSEDIRECTION = 5;
VOID = 6;
START = 7;
S1 = 8;
S2 = 9;
S3 = 10;
S4 = 11;
FOODWELL = 12;
END = 13;

CORRECT = 1;
INCORRECT = 2;

YBIN = 48; %the number of y bin
XBIN = 72; %the number of x bin
SCALE = 10; %pixel per bin.

areaN = 5; %area N : 1 = S12, 2 = S23, 3 = S34, 4 = intersection-foodwell, 5 = inbound

posFileName = 'VT1.nvt';
trialFileName = 'sessionSummary.csv';
posDataName = 'ParsedPosition';

posFile = [session_root '\' posFileName];
trialFile = [session_root '\ExtractedEvents\' trialFileName];

posDataDir = [session_root '\' posDataName];


%%-------------------------   variables initializing and explaning.    -----------------%%

[pData, ~] = simplePosPars(posFile, epochST, epochED); %position data
[tData, ~, ~] = trialPars(trialFile);
[trialN, ~] = size(tData{1,1});
[~, posN] = size(pData);

% voidN = sum(tData{1,VOID}(:,1));

%initialization
total_trial_number = trialN; % - voidN;
trial_void = logical(zeros(total_trial_number,1)); %it will be not used.
trial_correctness = zeros(total_trial_number,1);
trial_context = zeros(total_trial_number,1);
trial_side = zeros(total_trial_number,1);
trial_ambiguity = zeros(total_trial_number, 1);

a = zeros(posN, 1);
t = zeros(posN, 1);
x = zeros(posN, 1);
y = zeros(posN, 1);
trial = logical(zeros(posN, total_trial_number));
area = logical(zeros(posN, areaN));
correctness = logical(zeros(posN, 2));
cont = logical(zeros(posN, 6));
side = logical(zeros(posN, 2));
ambiguity = logical(zeros(posN, 3));
% void = logical(zeros(posN, 1));

PositionCountMap = zeros(YBIN, XBIN, total_trial_number);

%%---------------------       algorithm          ------------------------%%
%parsing trials. (total_trial_number, trial_corr, trial_context, trial_void)

trial_index = tData{1,INDEX};
trial_correctness = tData{1,CORRECTNESS};
trial_context = tData{1,CONTEXT};
trial_side = tData{1,RESPONSEDIRECTION};
trial_ambiguity = tData{1,AMBIGUITY};
 trial_void = tData{1,VOID};


 
%%
% tIter = 0; %real trial iteration. at first it exist because of void trial. 
%            %But trialPars program is revised. So there is no meaning of
%            %this variable.
% 
% for tempIter = 1:trialN
%      
%      %trial void
%      if tData{1,VOID}(tempIter,1) == 2
%          continue;
%      else
%          tIter = tIter + 1;
%          trial_void(tIter,1) = false;
%      end
%      
%      %trial correctness
%      if strcmp(tData{1,CORRECTNESS}(tempIter,1), 'TRUE')
%          trial_correctness(tIter,1) = CORRECT;
%      else
%          trial_correctness(tIter,1) = INCORRECT;
%      end
%      
%      contNum = tData{1,CONTEXT}(tempIter,1);
%      
%      %trial context
%      if contNum == 1
%          trial_context(tIter,1) = 1;
%      elseif contNum == 2
%          trial_context(tIter,1) = 2;
%      elseif contNum == 3
%          trial_context(tIter,1) = 3;
%      elseif contNum == 4
%          trial_context(tIter,1) = 4;     
%      else
%          error([sessionName ': this is not valid']);
%      end
%      
%  end
% 
%  total_trial_number = tIter;
 


%% parsing positions.(a, t, x, y ,corr, obj, side, trial, area)

cTrial = 1; %current trial.
cArea = 1; %current area.
validFlag = 0;

for pIter = 1:posN
    
    %head direction, time stamp, x y position.
    a(pIter,1) = pData(1,pIter).dir; %we will not use this variable)
    t(pIter,1) = pData(1,pIter).timeStamp / 1000000; %change the time unit
    x(pIter,1) = pData(1,pIter).xPos;
    y(pIter,1) = pData(1,pIter).yPos;
    
    %end of trial
    if cTrial > total_trial_number
        continue;
    end
    
    if validFlag == 0 && t(pIter,1) >= tData{1,START}(cTrial,1)     % new trial start
        msg = sprintf('%d %d %f %f\n', cTrial, pIter, t(pIter,1), tData{1,START}(cTrial,1));
%         disp(msg);
        validFlag = 1;
    end
    
    if validFlag == 1 && t(pIter,1) >= tData{1,END}(cTrial,1)       % trial end
        validFlag = 0;
        cArea = 1;
        cTrial = cTrial + 1;
    end
    
    if validFlag == 0
        continue;
    end
    
    %determine trial
    trial(pIter, cTrial) = true;
    
    %area determine
    if t(pIter,1) <= tData{1,S2}(cTrial, 1)
        cArea = 1;
    elseif t(pIter,1) <= tData{1,S3}(cTrial, 1)
        cArea = 2;
    elseif t(pIter,1) <= tData{1,S4}(cTrial, 1)
        cArea = 3;
    elseif t(pIter,1) <= tData{1,FOODWELL}(cTrial, 1)
        cArea = 4;
    else
        cArea = 5;
    end
    
    area(pIter, cArea) = true;
        
    correctness(pIter, trial_correctness(cTrial, 1)) = true;
    
    cont(pIter, trial_context(cTrial, 1)) = true;
    
    side(pIter, trial_side(cTrial, 1)) = true;
    
    ambiguity(pIter, trial_ambiguity(cTrial, 1)) = true;
    
%     void(pIter, trial_void(cTrial, 1)) = true;    
end


%make PositionCountMap

for tIter = 1:total_trial_number
    
    cLogicalTrial = logical(trial(:,tIter));
    if sum(cLogicalTrial)
        map = Pos2Map(x(cLogicalTrial), y(cLogicalTrial), XBIN, YBIN, SCALE);
        PositionCountMap(1:YBIN, 1:XBIN, tIter) = map;
    end
    
end

%% save parsed position results

save([session_root '\ParsedPosition.mat'], 't', 'x', 'y', 'a', 'cont', ...
    'correctness', 'trial', 'side', 'ambiguity', 'area', 'PositionCountMap', 'trial_correctness', 'trial_side', 'trial_context', 'trial_ambiguity', 'total_trial_number','-append'); %void trials are deleted.

if ~exist(posDataDir, 'dir')
    mkdir(posDataDir);
end

cd(posDataDir);

% display position for each trial (outbound + inbound)
tIter = 1;
while tIter <= total_trial_number
    h = figure;
    sheet_position = [2 2 25 25];
    set(gcf, 'units', 'centimeter', 'position', sheet_position, 'color', [1 1 1]);
    for iter = 1:16
        subplot(4,4,iter);
        plot(x(trial(:,tIter)), y(trial(:,tIter)), '.', 'MarkerSize', 3);
        set(gca, 'YDir', 'rev', 'XLim', [250 450], 'YLim', [50 500]);
        
        tIter = tIter + 1;
        
        if tIter > total_trial_number
            break;
        end
    end
    
    saveImage(h, ['total position from trial ' num2str(tIter) '.jpg'], 'centimeter', sheet_position);    
end

% display position for each trial (outbound only)
tIter = 1;
while tIter <= total_trial_number
    h = figure;
    sheet_position = [2 2 25 25];
    set(gcf, 'units', 'centimeter', 'position', sheet_position, 'color', [1 1 1]);
    for iter = 1:16
        h2 = subplot(4,4,iter); 
        plot(x(trial(:,tIter) & ~area(:,5)), y(trial(:,tIter) & ~area(:,5)), '.', 'MarkerSize', 3);
                
        set(gca, 'YDir', 'rev', 'XLim', [250 450], 'YLim', [50 500]);
        text(260,450,num2str(trial_index(tIter)),'Fontsize',12,'Parent',h2);

        tIter = tIter + 1;
        
        if tIter > total_trial_number
            break;
        end
    end
    
    saveImage(h, ['outbound position from trial ' num2str(tIter) '.jpg'], 'centimeter', sheet_position);    
end
