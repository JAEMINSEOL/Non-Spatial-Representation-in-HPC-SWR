
function createParsedSpike(mother_root, thisCLUSTER)

findHYPHEN = find(thisCLUSTER == '-');

thisRID = (thisCLUSTER(1, 1:findHYPHEN(1) - 1));
thisSID = (thisCLUSTER(1, findHYPHEN(1) + 1:findHYPHEN(2) - 1));
thisTTID = num2str(str2double((thisCLUSTER(1, findHYPHEN(2) + 1:findHYPHEN(3) - 1))));
thisCLID = num2str(str2double((thisCLUSTER(1, findHYPHEN(3) + 1:end))));

%Load Epoch information
cd([mother_root '\rat' thisRID]);
% if str2double(thisRID) > 400
%     thisEPOCH = ['behaviorEpoch_rat' thisRID '.xlsx'];
%     epochST = xlsread(thisEPOCH, ['A' thisSID ':' 'A' thisSID]);
%     epochED = xlsread(thisEPOCH, ['B' thisSID ':' 'B' thisSID]);
% else
    thisEPOCH = ['behaviorEpoch_rat' thisRID '.csv'];
    epochST = csvread(thisEPOCH, str2double(thisSID)-1, 0, [str2double(thisSID)-1, 0, str2double(thisSID)-1, 0]);
    epochED = csvread(thisEPOCH, str2double(thisSID)-1, 1, [str2double(thisSID)-1, 1, str2double(thisSID)-1, 1]);
% end

%Load cluster file
cd([mother_root '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID]);
if exist(['TT' thisTTID '_beh_SS_' thisCLID '.ntt'],'file')
    thisEpochCLTS = Nlx2MatSpike(['TT' thisTTID '_beh_SS_' thisCLID '.ntt'], [1 0 0 0 0], 0, 4, [epochST, epochED]);
    thisEpochCLTS = thisEpochCLTS ./ 1000000;
    spkN = length(thisEpochCLTS);
    area=[];
    %Load position parsing result
    load([mother_root '\rat' thisRID '\rat' thisRID '-' thisSID '\Behavior\ParsedPosition.mat']);
    areaN = 5;
    posN = length(x);
    if ~exist('ambiguity'), ambiguity = ones(posN,1); end
    if ~exist('cont'), cont = side; end
    
    a_spk = zeros(spkN, 1);
    t_spk = zeros(spkN, 1);
    x_spk = zeros(spkN, 1);
    y_spk = zeros(spkN, 1);
    trial_spk = logical(zeros(spkN, size(trial,2)));
    area_spk = logical(zeros(spkN, areaN));
    correctness_spk = logical(zeros(spkN, 2));
    cont_spk = logical(zeros(spkN, 6));
    side_spk = logical(zeros(spkN, 2));
    ambiguity_spk = logical(zeros(spkN, 3));
    %     void_spk = logical(zeros(spkN, 1));
    
    clRUN = 1;
    posRUN = 1;
    while posRUN <= posN
        if t(posRUN,1) >= thisEpochCLTS(clRUN) && posRUN == 1
            a_spk(clRUN) = a(posRUN);
            %         t_spk(clRUN) = t(posRUN);
            t_spk(clRUN) = thisEpochCLTS(clRUN);
            x_spk(clRUN) = x(posRUN);
            y_spk(clRUN) = y(posRUN);
            trial_spk(clRUN,:) = trial(posRUN,:);
            area_spk(clRUN,:) = area(posRUN,:);
            correctness_spk(clRUN,:) = correctness(posRUN,:);
            cont_spk(clRUN,:) = cont(posRUN,:);
            side_spk(clRUN,:) = side(posRUN,:);
            ambiguity_spk(clRUN,:) = ambiguity(posRUN,:);
            %         void_spk(clRUN) = void(posRUN);
            clRUN = clRUN + 1;
        elseif t(posRUN,1) >= thisEpochCLTS(clRUN)
            if t(posRUN,1)-thisEpochCLTS(clRUN) < thisEpochCLTS(clRUN)-t(posRUN-1,1)
                 a_spk(clRUN) = a(posRUN);
                %             t_spk(clRUN) = t(posRUN);
                t_spk(clRUN) = thisEpochCLTS(clRUN);
                x_spk(clRUN) = x(posRUN);
                y_spk(clRUN) = y(posRUN);
                trial_spk(clRUN,:) = trial(posRUN,:);
                area_spk(clRUN,:) = area(posRUN,:);
                correctness_spk(clRUN,:) = correctness(posRUN,:);
                cont_spk(clRUN,:) = cont(posRUN,:);
                side_spk(clRUN,:) = side(posRUN,:);
                ambiguity_spk(clRUN,:) = ambiguity(posRUN,:);
                %             void_spk(clRUN) = void(posRUN);
                posRUN = posRUN - 1; clRUN = clRUN + 1;
            else
                a_spk(clRUN) = a(posRUN-1);
                %             t_spk(clRUN) = t(posRUN-1);
                t_spk(clRUN) = thisEpochCLTS(clRUN);
                x_spk(clRUN) = x(posRUN-1);
                y_spk(clRUN) = y(posRUN-1);
                trial_spk(clRUN,:) = trial(posRUN-1,:);
                area_spk(clRUN,:) = area(posRUN-1,:);
                correctness_spk(clRUN,:) = correctness(posRUN-1,:);
                cont_spk(clRUN,:) = cont(posRUN-1,:);
                side_spk(clRUN,:) = side(posRUN-1,:);
                ambiguity_spk(clRUN,:) = ambiguity(posRUN-1,:);
                %             void_spk(clRUN) = void(posRUN-1);
                posRUN = posRUN - 1; clRUN = clRUN + 1;
            end
        end
        
        if clRUN > length(thisEpochCLTS)
            break;
        end
        
        posRUN = posRUN + 1;
    end
    
    if clRUN <= length(thisEpochCLTS)   % spikes after last position sampling.
        posRUN = posRUN - 1;
        
        for iter = clRUN : length(thisEpochCLTS)
            a_spk(iter) = a(posRUN);
            %         t_spk(clRUN) = t(posRUN);
            t_spk(iter) = thisEpochCLTS(iter);
            x_spk(iter) = x(posRUN);
            y_spk(iter) = y(posRUN);
            trial_spk(iter,:) = trial(posRUN,:);
            area_spk(iter,:) = area(posRUN,:);
            correctness_spk(iter,:) = correctness(posRUN,:);
            cont_spk(iter,:) = cont(posRUN,:);
            side_spk(iter,:) = side(posRUN,:);
            ambiguity_spk(iter,:) = ambiguity(posRUN,:);
            %         void_spk(clRUN) = void(posRUN);
            
        end
    end
    if strcmp(thisSID,'08')
     1;
    end
    save([mother_root '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' thisCLID '.mat'], 't_spk', 'x_spk', 'y_spk', 'a_spk', 'cont_spk', ...
        'correctness_spk', 'trial_spk', 'side_spk', 'ambiguity_spk', 'area_spk');
    
    figure;
     scatter(y,x,10,[.8 .8 .8],'filled')
axis off
hold on
scatter(y_spk,x_spk,10,'r','filled')
title(thisCLUSTER)
    
end % if exist
end