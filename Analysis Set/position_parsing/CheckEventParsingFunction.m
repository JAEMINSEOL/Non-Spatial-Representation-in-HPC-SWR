function CheckEventParsingFunction(session_root, thisRID, thisSID)

cd(session_root);
fid2 = fopen('behaviorData.csv','r');

cd([session_root,'\ExtractedEvents']);
fid1 = fopen('sessionSummary.csv','r');

max1 = 13; % number of row (sessionSummary)
max2 = 12; % number of row (behaviorData)

% index 1 for sessionSummary, index 2 for behaviorData

correctness1 = 2;
context1 = 3;
ambiguity1 = 4;
side1 = 5;
void1 = 6;
sensor1_1 = 8;
sensor2_1 = 9;
sensor3_1 = 10;
sensor4_1 = 11;
touch1 = 12;

start1 = 7; % for calculating latency
trial1 = 1; % for indicating

correctness2 = 4;
context2 = 2;
ambiguity2 = 12;
side2 = 3;
void2 = 10;
sensor1_2 = 9;
sensor2_2 = 8;
sensor3_2 = 7;
sensor4_2 = 6;
touch2 = 5;

% check session type
session_type = get_sessionType(thisRID, thisSID);
% 


% skip the headline
tline1 = fgetl(fid1);
tline2 = fgetl(fid2);

tline1 = fgetl(fid1);
tline2 = fgetl(fid2);

token1 = cell(max1,1);
token2 = cell(max2,1);

while ischar(tline1) && ischar(tline2)
    
    for iter = 1:max1
        [token1{iter}, tline1] = strtok(tline1,', ');
    end
    
    for iter = 1:max2
        [token2{iter}, tline2] = strtok(tline2,', ');
    end

    % change string to double
    token1{sensor1_1} = str2double(token1{sensor1_1});
    token1{sensor2_1} = str2double(token1{sensor2_1});    
    token1{sensor3_1} = str2double(token1{sensor3_1});    
    token1{sensor4_1} = str2double(token1{sensor4_1});    
    token1{touch1} = str2double(token1{touch1});    
    token1{start1} = str2double(token1{start1});
    
    token2{sensor1_2} = str2double(token2{sensor1_2});    
    token2{sensor2_2} = str2double(token2{sensor2_2});
    token2{sensor3_2} = str2double(token2{sensor3_2});
    token2{sensor4_2} = str2double(token2{sensor4_2});
    token2{touch2} = str2double(token2{touch2});
       
    % calculate the latencies
    token1{sensor1_1} = token1{sensor1_1} - token1{start1};
    token1{sensor2_1} = token1{sensor2_1} - token1{start1};
    token1{sensor3_1} = token1{sensor3_1} - token1{start1};
    token1{sensor4_1} = token1{sensor4_1} - token1{start1};
    token1{touch1} = token1{touch1} - token1{start1};
    
    % round off for time variables
%     token1{sensor1_1} = round(token1{sensor1_1}*10)/10;
%     token1{sensor2_1} = round(token1{sensor2_1}*10)/10;
%     token1{sensor3_1} = round(token1{sensor3_1}*10)/10;
%     token1{sensor4_1} = round(token1{sensor4_1}*10)/10;
%     token1{touch1} = round(token1{touch1}*10)/10;
%     token1{start1} = round(token1{start1}*10)/10;
%     
%     token2{sensor1_2} = round(token2{sensor1_2}*10)/10;    
%     token2{sensor2_2} = round(token2{sensor2_2}*10)/10;    
%     token2{sensor3_2} = round(token2{sensor3_2}*10)/10;    
%     token2{sensor4_2} = round(token2{sensor4_2}*10)/10;    
%     token2{touch2} = round(token2{touch2}*10)/10;
    

    % change characters to numbers
    switch token2{correctness2}
        case 'CORRECT'
            token2{correctness2} = '1';
        case 'WRONG'
            token2{correctness2} = '2';
    end
    
    switch token2{context2}
        case 'Zebra'
            token2{context2} = '1';
        case 'Pebbles'
            token2{context2} = '2';
        case 'Bamboo'
            token2{context2} = '3';
        case 'Mountain'
            token2{context2} = '4';
        case 'Peacock'
            token2{context2} = '5';
        case 'Palmtree'
            token2{context2} = '6';
        case 'Gray'
            token2{context2} = '7';
    end
    
    if session_type ~= 6
        
    switch token2{ambiguity2}
        case 'NORMAL'
            token2{ambiguity2} = '1';
        case 'AMB1'
            token2{ambiguity2} = '2';
        case 'AMB2'
            token2{ambiguity2} = '3';
        case 'SHORT'
            token2{ambiguity2} = '2';
        case 'MIDDLE'
            token2{ambiguity2} = '3';
    end
    
    end
    
    switch token2{side2}
        case 'left'
            token2{side2} = '1';
        case 'right'
            token2{side2} = '2';
    end
    
    switch token2{void2}
        case 'NO'
            token2{void2} = '1';
        case 'YES'
            token2{void2} = '2';
    end
    
    % compare the data
        
    if (token1{correctness1} ~= token2{correctness2}) && (token1{void1} == '1')
        disp(['mismatch occurs in ',session_root]);
        msg = sprintf(['trial = ',token1{trial1},' correctness mismatch! correctness1 = ',token1{correctness1},' correctness2 = ', token2{correctness2}]);  disp(msg);
    end
    if (token1{context1} ~= token2{context2}) && (token1{void1} == '1')
        disp(['mismatch occurs in ',session_root]);
        msg = sprintf(['trial = ',token1{trial1},' context mismatch! context1 = ',token1{context1},' context2 = ', token2{context2}]);  disp(msg);
    end
    
    if sum(session_type == [2 3 4])
        if (token1{ambiguity1} ~= token2{ambiguity2}) && (token1{void1} == '1')
            disp(['mismatch occurs in ',session_root]);
            msg = sprintf(['trial = ',token1{trial1},' ambiguity mismatch! ambiguity1 = ',token1{ambiguity1},' ambiguity2 = ', token2{ambiguity2}]);  disp(msg);
        end
    end
    
    if (token1{side1} ~= token2{side2}) && (token1{void1} == '1')
        disp(['mismatch occurs in ',session_root]);
        msg = sprintf(['trial = ',token1{trial1},' side mismatch! side1 = ',token1{side1},' side2 = ', token2{side2}]);  disp(msg);
    end
    if token1{void1} ~= token2{void2}
        disp(['mismatch occurs in ',session_root]);
        msg = sprintf(['trial = ',token1{trial1},' void mismatch! void1 = ',token1{void1},' void2 = ', token2{void2}]);  disp(msg);
    end
    if (abs(token1{sensor1_1} - token2{sensor1_2}) > 0.01) && (token1{void1} == '1')
        disp(['mismatch occurs in ',session_root]);
        msg = sprintf(['trial = ',token1{trial1},' sensor1 mismatch! sensor1_1 = %f, sensor1_2 = %f'],token1{sensor1_1}, token2{sensor1_2});  disp(msg);
    end
    if (abs(token1{sensor2_1} - token2{sensor2_2}) > 0.01) && (token1{void1} == '1')
        disp(['mismatch occurs in ',session_root]);
        msg = sprintf(['trial = ',token1{trial1},' sensor2 mismatch! sensor2_1 = %f, sensor2_2 = %f'],token1{sensor2_1}, token2{sensor2_2});  disp(msg);
    end
    if (abs(token1{sensor3_1} - token2{sensor3_2}) > 0.01) && (token1{void1} == '1')
        disp(['mismatch occurs in ',session_root]);
        msg = sprintf(['trial = ',token1{trial1},' sensor3 mismatch! sensor3_1 = %f, sensor3_2 = %f'],token1{sensor3_1}, token2{sensor3_2});  disp(msg);
    end
    if (abs(token1{sensor4_1} - token2{sensor4_2}) > 0.01) && (token1{void1} == '1')
        disp(['mismatch occurs in ',session_root]);
        msg = sprintf(['trial = ',token1{trial1},' sensor4 mismatch! sensor4_1 = %f, sensor4_2 = %f'],token1{sensor4_1}, token2{sensor4_2});  disp(msg);
    end
    if (abs(token1{touch1} - token2{touch2}) > 0.01) && (token1{void1} == '1')
        disp(['mismatch occurs in ',session_root]);
        msg = sprintf(['trial = ',token1{trial1},' touch mismatch! touch1 = %f, touch2 = %f'],token1{touch1}, token2{touch2});  disp(msg);
    end
   
    tline1 = fgetl(fid1);
    tline2 = fgetl(fid2);
    
end

if ischar(tline1) || ischar(tline2)
    disp(['mismatch occurs in ',session_root,'. Line number mismatch!']);
end

disp('checking finished');

end