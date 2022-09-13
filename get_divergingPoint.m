
% 2015-Jun-18.
% 2022-Sep-13 (SJM). Add "return 0" in line 32

function diverging_point = get_divergingPoint(mother_root, thisRID, thisSID, opt)

diverging_point = [];

if nargin == 3 || (nargin == 4 && strcmp(opt,'outbound'))    
    fid = fopen([mother_root '\diverging_points.csv'], 'r');    
elseif nargin == 4 && strcmp(opt,'inbound')    
    fid = fopen([mother_root '\diverging_points_inbound.csv'], 'r');    
end

fline = fgetl(fid);  % read header

fline = fgetl(fid);
while ischar(fline)
    [session_temp, diverging_point_temp] = strtok(fline, ',');
    
    if strcmp(session_temp, [thisRID '-' thisSID]);
        diverging_point = str2num(diverging_point_temp);
        break;
    end
    fline = fgetl(fid);
end

fclose(fid);

if size(diverging_point) == 0
    disp('This session has no diverging point!');
    diverging_point = 0;
end

end