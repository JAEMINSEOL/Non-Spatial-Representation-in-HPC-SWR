
function EventDecoderFunction(session_root, end2Flag, track_type)

cd(session_root);
evtfile_listing = dir('Events.txt');

fid = fopen('Events.txt','r');
if ~exist('ExtractedEvents','dir')
    mkdir('ExtractedEvents');
end
cd('ExtractedEvents');
%     listing=dir('trial*.mat');[r_trial,c_trial]=size(listing);
% %     if r_trial, return; end

% prepare saving xls file
xlsName = 'sessionSummary.csv';
xlsHEADER = 'Trial#, correctness, context, ambiguity, response side, void, start, sensor1, sensor2, sensor3, sensor4, food sensor, end\n';
xlsID = fopen(xlsName,'W');

fprintf(xlsID, xlsHEADER);


n=evtfile_listing(1).bytes;
n2=0;
time=0;go=true ;
directionFlag = 0; % 1 for go, 0 for return

% there's no header line in our events files, so just skip
% tline = fgetl(fid); %header
% n2=n2+length(tline);waitbar(n2/n);
% tline = fgetl(fid); %line jump
% n2=n2+length(tline);waitbar(n2/n);

h=waitbar(0);
waitbar(n2/n);

t_line=1;

% reset all variables
correctness=[]; context=[]; ambiguity=1; side=[]; void=1;
sensor1=[]; sensor2=[]; sensor3=[]; sensor4=[];
start=[]; touch=[]; stop=[];
zone0=[]; zone1=[]; zone2=[]; zone3=[]; zone4=[];

trial_number=0;tline = fgetl(fid);n2=n2+length(tline);waitbar(n2/n);go=1;

start_number=0;end_number=0;


while ischar(tline)     % find 'behavior recording start' event
    
    cor = 1;
    
    for cor = 1:18
        [token, tline] = strtok(tline,',');
    end
    
    if strfind(token,'beh')
        break;
    end
    
    tline = fgetl(fid);
    n2 = n2+length(tline);
    waitbar(n2/n);
    t_line = t_line+1;
end

tline = fgetl(fid);
n2 = n2+length(tline);
waitbar(n2/n);
t_line = t_line+1;

p_array =[];
index = [1:8];
while ischar(tline)
    
    for col = 1:4
        [token, tline] = strtok(tline,',');
    end
    time = str2double(token)/1000000;
    
    for col = 5:18
        [token, tline] = strtok(tline,',');
    end
    
    if strfind(token,'Zoned Video:')
        if directionFlag
            a = token(19);  % zone number
            p = token(22);  % p = 'n' : zone entered, p = 'x' : zone exited
            
            if p == 'n'
                p = 1;
            elseif p == 'x'
                p = 2;
            end
            
            if a == '0'
                zone0(p,end+1) = time;
            end
            if a == '1'
                zone1(p,end+1) = time;
            end
            if a == '2'
                zone2(p,end+1) = time;
            end
            if a == '3'
                zone3(p,end+1) = time;
            end
            if a == '4'
                zone4(p,end+1) = time;
            end
        end
    end
    
    if strfind(token,'TTL Input on AcqSystem1_0 board 0 port')
        port = token(41);
        a = token(54:55);
        p = dec2binvec( hex2dec(a), 8);
        
        p_array(end+1,index) = p;
        
        %         p_diff = xor(p_old,p);
        %         p_old = p;
        
        if port == '1' || port == '3'        % sensor related signal
            if p(5) == 0
                sensor1(end+1) = time;
            end
            if p(6) && ~isempty(sensor1)
                sensor2(end+1) = time;
            end
            if p(7) && ~isempty(sensor2)
                sensor3(end+1) = time;
            end
            if p(8) %&& ~isempty(sensor3)
                sensor4(end+1) = time;
            end
        elseif port == '0' || port == '2'
            if p(2) && p(3)     % photosensor detecting time at task 1 (mask signal)
                if isempty(touch)
                    touch = time;
                    %                 msg=sprintf('touch t=%1.2f ',touch); disp(msg);
                    directionFlag = 0;
                end
            end
            
            if p(2) == 1 && (p(3) == 0 && p(4) == 0)          % start mask signal detected,
                
                % check whether this is real start signal or not.
                % real one is followed by commit signal.
                tline = fgetl(fid);
                n2 = n2+length(tline);
                waitbar(n2/n);
                t_line = t_line+1;
                for col = 1:18
                    [token, tline] = strtok(tline,',');
                end
                
                if token(2) == 'T' && (token(41) == '0' || token(41) == '2')
                    
                    if ~(start_number == end_number)          % signal non-match (end signal missing)
                        disp('trial data(end signal) are missing');
                        min_=min([start_number end_number]);
                        max_=max([start_number end_number]);
                        %sensor1=[]; sensor2=[];sensor3=[]; sensor4=[];obj_touch=[];touch=[];trial_void=0;object=[];correctness=[];side=[];start=[];stop=[];
                        sensor1=[]; sensor2=[];sensor3=[]; sensor4=[]; context=[]; side=[]; correctness=[]; touch=[]; void=0; stop=[];     % only keep start value
                        zone0=[]; zone1=[]; zone2=[]; zone3=[]; zone4=[];
                        
                        if ~min_, min_=1; end
                        
                        for ii=min_:(max_ -1)
                            
                            nb_digit= numel(num2str(ii));
                            
                            switch nb_digit
                                
                                case 1
                                    f= sprintf('trial000%d.mat',ii);
                                case 2
                                    f= sprintf('trial00%d.mat',ii);
                                case 3
                                    f = sprintf('trial0%d.mat',ii);
                                otherwise
                            end
                            
                            void=1;
                            %save(f,'sensor1','sensor2','sensor3','sensor4','obj_touch','touch','trial_void','object','correctness','side','start','stop');
                            save(f,'sensor1','sensor2','sensor3','sensor4','zone0','zone1','zone2','zone3','zone4','start','context','side','correctness','touch','void','stop');
                        end
                        start_number = min_;
                        end_number=min_;
                        
                        trial_number = trial_number - 1;
                    end
                    
                    % reset all variables
                    correctness=[]; context=[]; ambiguity=1; side=[]; void=1;
                    sensor1=[]; sensor2=[]; sensor3=[]; sensor4=[];
                    start=[]; touch=[]; stop=[];
                    zone0=[]; zone1=[]; zone2=[]; zone3=[]; zone4=[];
                    
                    trial_number = trial_number+1;
                    start_number = start_number+1;
                    start = time;
                    %                 msg = sprintf('start #%d, t=%1.2f',trial_number,time); disp(msg);
                    directionFlag = 1;
                end
            end
            
            if (p(3) == 1) && (p(2) == 0 && p(4) == 0)          % end 1 signal detected
                end_number = end_number+1;
                stop = time;
                %                 msg = sprintf('stop t=%1.2f',time); disp(msg);
                
                % correctness
                % correct = 1, incorrect = 2;
                if p(5)
                    correctness = 1;
                else
                    correctness = 2;
                end
                %                 msg = sprintf('correctness = %d',correctness); disp(msg);
                
                % void
                % not void trial = 1, void trial = 2;
                if p(6)
                    void = 2;
                    %                     msg = sprintf('VOID');disp(msg);
                else
                    void = 1;
                end
                
                % context
                if p(7) == 0 && p(8) == 0   % zebra
                    context = 1;
                elseif p(7) == 0 && p(8) == 1   % pebble
                    context = 2;
                elseif p(7) == 1 && p(8) == 0   % bamboo
                    context = 3;
                    %                     context = 5;
                elseif p(7) == 1 && p(8) == 1   % mountain
                    context = 4;
                    %                     context = 6;
                end
                
                
                %                 switch context
                %                     case 1
                %                         msg = sprintf('context : zebra'); disp(msg);
                %                     case 2
                %                         msg = sprintf('context : pebble'); disp(msg);
                %                     case 3
                %                         msg = sprintf('context : bamboo'); disp(msg);
                %                     case 4
                %                         msg = sprintf('context : mountain'); disp(msg);
                %                 end
                
                % choice side
                % side : left=1, right=2;
                if mod(context, 2) == 1 && p(5) == 1
                    side = 1;
                elseif mod(context, 2) == 1 && p(5) == 0
                    side = 2;
                elseif mod(context, 2) == 0 && p(5) == 1
                    side = 2;
                elseif mod(context, 2) == 0 && p(5) == 0
                    side = 1;
                end
                
                
                % check whether there is end2 signal.
                % end2 signal stands for 3rd pair of context(end2Flag = 1) or
                % ambiguity/trace(end2Flag = 2).
                if end2Flag ~= 0
                    
                    tline = fgetl(fid);
                    n2 = n2+length(tline);
                    waitbar(n2/n);
                    t_line = t_line+1;
                    for col = 1:18
                        [token, tline] = strtok(tline,',');
                    end
                    
                    port = token(41);
                    a = token(54:55);
                    p = dec2binvec( hex2dec(a), 8);
                    
                    if end2Flag == 1
                        if port == '0' && p(4) == 1
                            if p(8)
                                context = 5;
                            elseif p(7)
                                context = 6;
                                side = 3-side;
                            end
                        end
                    elseif end2Flag == 2
                        if p(6) == 1 || p(8) == 1,  ambiguity = 2;
                        elseif p(5) == 1 || p(7) == 1,  ambiguity = 3;
                        end
                    end
                end
                
                if (start_number == end_number)
                    if isempty(touch) % if it is skip trial (so no signal form foodwell sensor), touch time = trial end time
                        touch = stop;
                    end
                    
                    if track_type == 1, sensor4(1) = 0; end
                    
                    f=sprintf('trial%d.mat',trial_number);
                    save(f,'start_number','end_number','sensor1','sensor2','sensor3','sensor4','zone0','zone1','zone2','zone3','zone4','start','context','ambiguity','side','correctness','touch','void','stop');
                    
                    % save trial summary xls file
                    if isempty(sensor3)
                        sensor3 = -1;
                    end
                    if isempty(sensor4)
                        sensor4 = -1;
                    end
                    result = [trial_number, correctness, context, ambiguity, side, void, start, sensor1(1), sensor2(1), sensor3(1), sensor4(1), touch, stop];
                    fprintf(xlsID, '%d, %d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f\n', result);
                    
                else        % signal non-match (start signal missing)
                    disp('trial data(start signal) are missing');
                    min_=min([start_number end_number]);
                    max_=max([start_number end_number]);
                    sensor1=[]; sensor2=[]; sensor3=[];sensor4=[]; start=[]; context=[]; side=[]; correctness=[]; touch=[]; void=0;        % keep stop value
                    zone0=[]; zone1=[]; zone2=[]; zone3=[]; zone4=[];
                    if ~min_, min_=1; end
                    for ii=min_:(max_ - 1)
                        
                        nb_digit= numel(num2str(ii));
                        
                        switch nb_digit
                            
                            case 1
                                f= sprintf('trial000%d.mat',ii);
                            case 2
                                f= sprintf('trial00%d.mat',ii);
                            case 3
                                f = sprintf('trial0%d.mat',ii);
                            otherwise
                        end
                        
                        void=1;
                        save(f,'sensor1','sensor2','sensor3','sensor4','zone0','zone1','zone2','zone3','zone4','start','context','ambiguity', 'side','correctness','touch','void','stop');
                    end
                    end_number=min_;
                    start_number = min_;
                end
                
            elseif p(5)     % end signal for task 2
            end
        end
        
    elseif strfind(token,'Digital Lynx Parallel Input Port TTL')
        % 0x8211 = sensor 1(stbox) 0x8201-0x9201
        % sensor 2 이후 : 0x1201-0x5201
        port = token(41);
        a = token(42:43);

        
        %         p_diff = xor(p_old,p);
        %         p_old = p;
        
        if port == '1' || port == '3'        % sensor related signal
            if p(5) == 0
                sensor1(end+1) = time;
            end
            if p(6) && ~isempty(sensor1)
                sensor2(end+1) = time;
            end
            if p(7) && ~isempty(sensor2)
                sensor3(end+1) = time;
            end
            if p(8) %&& ~isempty(sensor3)
                sensor4(end+1) = time;
            end
        elseif port == '0' || port == '2'
            if p(2) && p(3)     % photosensor detecting time at task 1 (mask signal)
                if isempty(touch)
                    touch = time;
                    %                 msg=sprintf('touch t=%1.2f ',touch); disp(msg);
                    directionFlag = 0;
                end
            end
            
            if p(2) == 1 && (p(3) == 0 && p(4) == 0)          % start mask signal detected,
                
                % check whether this is real start signal or not.
                % real one is followed by commit signal.
                tline = fgetl(fid);
                n2 = n2+length(tline);
                waitbar(n2/n);
                t_line = t_line+1;
                for col = 1:18
                    [token, tline] = strtok(tline,',');
                end
                
                if token(2) == 'T' && (token(41) == '0' || token(41) == '2')
                    
                    if ~(start_number == end_number)          % signal non-match (end signal missing)
                        disp('trial data(end signal) are missing');
                        min_=min([start_number end_number]);
                        max_=max([start_number end_number]);
                        %sensor1=[]; sensor2=[];sensor3=[]; sensor4=[];obj_touch=[];touch=[];trial_void=0;object=[];correctness=[];side=[];start=[];stop=[];
                        sensor1=[]; sensor2=[];sensor3=[]; sensor4=[]; context=[]; side=[]; correctness=[]; touch=[]; void=0; stop=[];     % only keep start value
                        zone0=[]; zone1=[]; zone2=[]; zone3=[]; zone4=[];
                        
                        if ~min_, min_=1; end
                        
                        for ii=min_:(max_ -1)
                            
                            nb_digit= numel(num2str(ii));
                            
                            switch nb_digit
                                
                                case 1
                                    f= sprintf('trial000%d.mat',ii);
                                case 2
                                    f= sprintf('trial00%d.mat',ii);
                                case 3
                                    f = sprintf('trial0%d.mat',ii);
                                otherwise
                            end
                            
                            void=1;
                            %save(f,'sensor1','sensor2','sensor3','sensor4','obj_touch','touch','trial_void','object','correctness','side','start','stop');
                            save(f,'sensor1','sensor2','sensor3','sensor4','zone0','zone1','zone2','zone3','zone4','start','context','side','correctness','touch','void','stop');
                        end
                        start_number = min_;
                        end_number=min_;
                        
                        trial_number = trial_number - 1;
                    end
                    
                    % reset all variables
                    correctness=[]; context=[]; ambiguity=1; side=[]; void=1;
                    sensor1=[]; sensor2=[]; sensor3=[]; sensor4=[];
                    start=[]; touch=[]; stop=[];
                    zone0=[]; zone1=[]; zone2=[]; zone3=[]; zone4=[];
                    
                    trial_number = trial_number+1;
                    start_number = start_number+1;
                    start = time;
                    %                 msg = sprintf('start #%d, t=%1.2f',trial_number,time); disp(msg);
                    directionFlag = 1;
                end
            end
            
            if (p(3) == 1) && (p(2) == 0 && p(4) == 0)          % end 1 signal detected
                end_number = end_number+1;
                stop = time;
                %                 msg = sprintf('stop t=%1.2f',time); disp(msg);
                
                % correctness
                % correct = 1, incorrect = 2;
                if p(5)
                    correctness = 1;
                else
                    correctness = 2;
                end
                %                 msg = sprintf('correctness = %d',correctness); disp(msg);
                
                % void
                % not void trial = 1, void trial = 2;
                if p(6)
                    void = 2;
                    %                     msg = sprintf('VOID');disp(msg);
                else
                    void = 1;
                end
                
                % context
                if p(7) == 0 && p(8) == 0   % zebra
                    context = 1;
                elseif p(7) == 0 && p(8) == 1   % pebble
                    context = 2;
                elseif p(7) == 1 && p(8) == 0   % bamboo
                    context = 3;
                    %                     context = 5;
                elseif p(7) == 1 && p(8) == 1   % mountain
                    context = 4;
                    %                     context = 6;
                end
                
                
                %                 switch context
                %                     case 1
                %                         msg = sprintf('context : zebra'); disp(msg);
                %                     case 2
                %                         msg = sprintf('context : pebble'); disp(msg);
                %                     case 3
                %                         msg = sprintf('context : bamboo'); disp(msg);
                %                     case 4
                %                         msg = sprintf('context : mountain'); disp(msg);
                %                 end
                
                % choice side
                % side : left=1, right=2;
                if mod(context, 2) == 1 && p(5) == 1
                    side = 1;
                elseif mod(context, 2) == 1 && p(5) == 0
                    side = 2;
                elseif mod(context, 2) == 0 && p(5) == 1
                    side = 2;
                elseif mod(context, 2) == 0 && p(5) == 0
                    side = 1;
                end
                
                
                % check whether there is end2 signal.
                % end2 signal stands for 3rd pair of context(end2Flag = 1) or
                % ambiguity/trace(end2Flag = 2).
                if end2Flag ~= 0
                    
                    tline = fgetl(fid);
                    n2 = n2+length(tline);
                    waitbar(n2/n);
                    t_line = t_line+1;
                    for col = 1:18
                        [token, tline] = strtok(tline,',');
                    end
                    
                    port = token(41);
                    a = token(54:55);
                    p = dec2binvec( hex2dec(a), 8);
                    
                    if end2Flag == 1
                        if port == '0' && p(4) == 1
                            if p(8)
                                context = 5;
                            elseif p(7)
                                context = 6;
                                side = 3-side;
                            end
                        end
                    elseif end2Flag == 2
                        if p(6) == 1 || p(8) == 1,  ambiguity = 2;
                        elseif p(5) == 1 || p(7) == 1,  ambiguity = 3;
                        end
                    end
                end
                
                if (start_number == end_number)
                    if isempty(touch) % if it is skip trial (so no signal form foodwell sensor), touch time = trial end time
                        touch = stop;
                    end
                    
                    if track_type == 1, sensor4(1) = 0; end
                    
                    f=sprintf('trial%d.mat',trial_number);
                    save(f,'start_number','end_number','sensor1','sensor2','sensor3','sensor4','zone0','zone1','zone2','zone3','zone4','start','context','ambiguity','side','correctness','touch','void','stop');
                    
                    % save trial summary xls file
                    if isempty(sensor3)
                        sensor3 = -1;
                    end
                    if isempty(sensor4)
                        sensor4 = -1;
                    end
                    result = [trial_number, correctness, context, ambiguity, side, void, start, sensor1(1), sensor2(1), sensor3(1), sensor4(1), touch, stop];
                    fprintf(xlsID, '%d, %d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %f\n', result);
                    
                else        % signal non-match (start signal missing)
                    disp('trial data(start signal) are missing');
                    min_=min([start_number end_number]);
                    max_=max([start_number end_number]);
                    sensor1=[]; sensor2=[]; sensor3=[];sensor4=[]; start=[]; context=[]; side=[]; correctness=[]; touch=[]; void=0;        % keep stop value
                    zone0=[]; zone1=[]; zone2=[]; zone3=[]; zone4=[];
                    if ~min_, min_=1; end
                    for ii=min_:(max_ - 1)
                        
                        nb_digit= numel(num2str(ii));
                        
                        switch nb_digit
                            
                            case 1
                                f= sprintf('trial000%d.mat',ii);
                            case 2
                                f= sprintf('trial00%d.mat',ii);
                            case 3
                                f = sprintf('trial0%d.mat',ii);
                            otherwise
                        end
                        
                        void=1;
                        save(f,'sensor1','sensor2','sensor3','sensor4','zone0','zone1','zone2','zone3','zone4','start','context','ambiguity', 'side','correctness','touch','void','stop');
                    end
                    end_number=min_;
                    start_number = min_;
                end
                
            elseif p(5)     % end signal for task 2
            end
        end
        
    end
    
    tline = fgetl(fid);
    n2 = n2+length(tline);
    waitbar(n2/n);
    t_line = t_line+1;
end

fclose(xlsID);






