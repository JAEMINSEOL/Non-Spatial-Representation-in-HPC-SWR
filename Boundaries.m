
% 2015-Jun-19
% Boundary

% Boundaries from separated_condition (2015-Jan-20 ver.) %
% use switch-case and add default diverging_point value (2022-Sep-13 ver., SJM)

% 1 for short track, 2 for long track.
if strcmp(thisRID, '214') || strcmp(thisRID, '219'), track_type = 1;
elseif str2double(thisRID) < 400, track_type = 2;
elseif str2double(thisRID) > 400 && str2double(thisRID) < 500, track_type = 4;
elseif str2double(thisRID) > 500, track_type = 3;
end
%

switch str2double(thisRID)
    case {214 219} % LSM T-maze
        xEdge = [340 400 400 420 420 300 300 340 340];
        yEdge = [390 390 230 230 185 185 230 230 390];
        xEdge = xEdge*0.23;
        yEdge = yEdge*0.23;
        imROW = 100;
        imCOL = 130;
        
    case {232 234 295 313} % LSM T-maze
        xEdge = [325 395 395 475 475 245 245 325 325];
        yEdge = [480 480 160 160 70 70 160 160 480];
        xEdge = xEdge*0.23;
        yEdge = yEdge*0.23;
        imROW = 100;
        imCOL = 130;
        
    case {553 561 562} % LSM & LCH T-maze
        xEdge = [325 395 395 475 475 245 245 325 325];
        yEdge = [480 480 200 200 110 110 200 200 490];
        xEdge = xEdge*0.23;
        yEdge = yEdge*0.23;
        imROW = 100;
        imCOL = 130;
        
    case {415 425 427 454 471 487} % LSM & LCH T-maze
        xEdge = [325 395 395 475 475 245 245 325 325];
        yEdge = [480 480 190 190 100 100 190 190 480];
        xEdge = xEdge*0.23;
        yEdge = yEdge*0.23;
        imROW = 100;
        imCOL = 130;
        
    case {477 526 536 551 564 578 588} % Jhoseph 1D-VR rats
        xEdge = [0 2 2 0 ];
        yEdge = [300 300 0 0];
        imROW = 800;
        imCOL = 200;
        
        [y,x] = deal(x,y);
        [y_spk,x_spk] = deal(x_spk, y_spk);
        
        
    case {027 037 049 052 057 069 080} % Delcasso touchscreen rats
        xEdge = [300 380 380 430 430 250 250 300 300];
        yEdge = [400 400 200 200 150 150 200 200 400];
        xEdge = xEdge*0.23;
        yEdge = yEdge*0.23;
        imROW = 100;
        imCOL = 130;
    otherwise % rat >500
        xEdge = [325 395 395 475 475 245 245 325 325];
        yEdge = [480 480 200 200 110 110 200 200 490];
        xEdge = xEdge*0.23;
        yEdge = yEdge*0.23;
        imROW = 100;
        imCOL = 130;
        
end

if diverging_point==0, diverging_point = min(yEdge); end
stem_boundary = [xEdge(1), xEdge(1), xEdge(2), xEdge(2); yEdge(1), diverging_point, diverging_point, yEdge(1)];

% Corner_boundary;  % copy from separatedCondition (2014-Dec-23 ver.)

% Corner_boundary
% This specifies corner boundaries for each track sessions.
% Position traces of short track sessions were too different each other, so
% they needed their own corner boundary.
%Code by Hyunwoo Lee, 2014-Dec-25


if strcmp(thisRID, '214')
    if strcmp(thisSID, '05')
        corner_boundary_left = [340, 340, 375, 375; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [345, 345, 380, 380; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '06')
        corner_boundary_left = [340, 340, 375, 375; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [345, 345, 380, 380; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '07')
        corner_boundary_left = [340, 340, 380, 380; diverging_point, 190, 190, diverging_point];
        corner_boundary_right = [345, 345, 385, 385; diverging_point, 190, 190, diverging_point];
    elseif strcmp(thisSID, '08')
        corner_boundary_left = [345, 345, 375, 375; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [350, 350, 375, 375; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '12')
        corner_boundary_left = [335, 335, 365, 365; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [335, 335, 375, 375; diverging_point, 195, 195, diverging_point];
    else
        corner_boundary_left = [340, 340, 380, 380; diverging_point, 190, 190, diverging_point];
        corner_boundary_right = [350, 350, 390, 390; diverging_point, 190, 190, diverging_point];
    end
    
elseif strcmp(thisRID, '219')
    if strcmp(thisSID, '08')
        corner_boundary_left = [340, 340, 380, 380; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [350, 350, 380, 380; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '09')
        corner_boundary_left = [345, 345, 390, 390; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [345, 345, 380, 380; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '10')
        corner_boundary_left = [340, 340, 380, 380; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [345, 345, 385, 385; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '11')
        corner_boundary_left = [340, 340, 385, 385; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [345, 345, 385, 385; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '12')
        corner_boundary_left = [340, 340, 375, 375; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [345, 345, 385, 385; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '13')
        corner_boundary_left = [340, 340, 380, 380; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [350, 350, 385, 385; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '14')
        corner_boundary_left = [345, 345, 385, 385; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [350, 350, 390, 390; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '15')
        corner_boundary_left = [345, 345, 380, 380; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [355, 355, 385, 385; diverging_point, 195, 195, diverging_point];
    elseif strcmp(thisSID, '16') || strcmp(thisSID, '17')
        corner_boundary_left = [345, 345, 385, 385; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [350, 350, 385, 385; diverging_point, 195, 195, diverging_point];
    else
        corner_boundary_left = [340, 340, 380, 380; diverging_point, 195, 195, diverging_point];
        corner_boundary_right = [350, 350, 390, 390; diverging_point, 195, 195, diverging_point];
    end
    
elseif strcmp(thisRID, '232')
    corner_boundary_left = [335, 335, 390, 390; 70, diverging_point, diverging_point, 70];
    corner_boundary_right = [320, 320, 380, 380; 70, diverging_point, diverging_point, 70];
    
elseif strcmp(thisRID, '234')
    if strcmp(thisSID, '01') || strcmp(thisSID, '03')
        corner_boundary_left = [330, 330, 390, 390; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [320, 320, 380, 380; 70, diverging_point, diverging_point, 70];
    elseif strcmp(thisSID, '04')
        corner_boundary_left = [330, 330, 390, 390; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [320, 320, 390, 390; 70, diverging_point, diverging_point, 70];
    else
        corner_boundary_left = [335, 335, 390, 390; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [320, 320, 380, 380; 70, diverging_point, diverging_point, 70];
    end
    
elseif strcmp(thisRID, '285')
    if strcmp(thisSID, '03') || strcmp(thisSID, '04')
        corner_boundary_left = [335, 335, 390, 390; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [320, 320, 375, 375; 70, diverging_point, diverging_point, 70];
    else
        corner_boundary_left = [335, 335, 390, 390; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [320, 320, 380, 380; 70, diverging_point, diverging_point, 70];
    end
    
elseif strcmp(thisRID, '295')
    if strcmp(thisSID, '02') || strcmp(thisSID, '03') || strcmp(thisSID, '04') || strcmp(thisSID, '05')
        corner_boundary_left = [330, 330, 390, 390; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [320, 320, 380, 380; 70, diverging_point, diverging_point, 70];
    elseif strcmp(thisSID, '07') || strcmp(thisSID, '09')
        corner_boundary_left = [335, 335, 390, 390; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [320, 320, 370, 370; 70, diverging_point, diverging_point, 70];
    elseif strcmp(thisSID, '12')
        corner_boundary_left = [320, 320, 405, 405; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [300, 300, 385, 385; 70, diverging_point, diverging_point, 70];
    else
        corner_boundary_left = [335, 335, 390, 390; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [320, 320, 380, 380; 70, diverging_point, diverging_point, 70];
    end
    
elseif strcmp(thisRID, '415')
    if strcmp(thisSID, '14')
        corner_boundary_left = [325, 325, 430, 430; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [280, 280, 395, 395; 70, diverging_point, diverging_point, 70];
    elseif strcmp(thisSID, '16')
        corner_boundary_left = [325, 325, 430, 430; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [280, 280, 395, 395; 70, diverging_point, diverging_point, 70];
    else
        %          corner_boundary_left = [295, 295, 450, 450; 70, diverging_point, diverging_point, 70];
        %          corner_boundary_right = [260, 260, 420, 420; 70, diverging_point, diverging_point, 70];
        corner_boundary_left = [325, 325, 430, 430; 70, diverging_point, diverging_point, 70];
        corner_boundary_right = [280, 280, 395, 395; 70, diverging_point, diverging_point, 70];
    end
    
elseif strcmp(thisRID, '427')
    corner_boundary_left = [320, 320, 400, 400; 70, diverging_point, diverging_point, 70];
    corner_boundary_right = [300, 300, 380, 380; 70, diverging_point, diverging_point, 70];
    %     if strcmp(thisSID, '04')
    %         corner_boundary_left = [320, 320, 390, 390; 70, diverging_point, diverging_point, 70];
    %         corner_boundary_right = [320, 320, 380, 380; 70, diverging_point, diverging_point, 70];
    %     else
    %         corner_boundary_left = [335, 335, 390, 390; 70, diverging_point, diverging_point, 70];
    %         corner_boundary_right = [320, 320, 380, 380; 70, diverging_point, diverging_point, 70];
    %     end
    
elseif strcmp(thisRID, '561')
    corner_boundary_left = [330, 330, 400, 400; 70, diverging_point, diverging_point, 70];
    corner_boundary_right = [300, 300, 380, 380; 70, diverging_point, diverging_point, 70];
else

    corner_boundary_left = [xEdge(1), xEdge(1), xEdge(2), xEdge(2); min(yEdge), diverging_point, diverging_point, min(yEdge)]/0.23;
    corner_boundary_right = [xEdge(1), xEdge(1), xEdge(2), xEdge(2); min(yEdge), diverging_point, diverging_point, min(yEdge)]/0.23;
end

corner_boundary_left = corner_boundary_left*0.23;
corner_boundary_right = corner_boundary_right*0.23;