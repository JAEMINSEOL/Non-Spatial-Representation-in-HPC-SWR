

function [field_count, start_index, end_index, field_size, h] = field_boundary_function_jm(overall_map, threshold_ratio)

baseline_fr = 0;
    
% Find nan bins
nan_index(1) = min(find(~isnan(overall_map) == true)) - 1;
nan_index(2) = max(find(~isnan(overall_map) == true)) + 1;
overall_map(isnan(overall_map)) = [];

map_length = length(overall_map);
%

% %

overall_map = overall_map - baseline_fr;
overall_map(overall_map < 0) = 0;

% % Find peaks

% overall_map_smooth = smooth(overall_map);
overall_map_smooth = overall_map;
max_rates = max(overall_map_smooth);
if max_rates * 0.3 < 0.5, min_prominence = 0.5;
else, min_prominence = max_rates * 0.3;
end
[peak_rates, peak_locs] = findpeaks(overall_map_smooth, 'MinPeakProminence', min_prominence);


middle_peaknum = length(peak_locs);


% Check left end
temp = [];
temp(1) = 0;
temp(2 : length(overall_map_smooth) + 1) = overall_map_smooth;

[temp_rates, temp_locs, temp_width, temp_proms] = findpeaks(temp, 'MinPeakProminence', min_prominence);

if length(temp_locs) > middle_peaknum
    peak_rates(end+1) = temp_rates(1);
    peak_locs(end+1) = temp_locs(1) - 1; if peak_locs(end) < 1, peak_locs(end) = 1; end
end
%

% Check right end
temp = [];
temp(1 : length(overall_map_smooth)) = overall_map_smooth;
temp(end+1) = 0;

[temp_rates, temp_locs, temp_width, temp_proms] = findpeaks(temp, 'MinPeakProminence', min_prominence);

if length(temp_locs) > middle_peaknum
    peak_rates(end+1) = temp_rates(end);
    peak_locs(end+1) = temp_locs(end);
end
% %



% % Set field boundary

for iter = 1 : length(peak_locs)
    
    threshold_rate = peak_rates(iter) * threshold_ratio;
    
    % Find start point of the field
    flag = 1;
    current_index = peak_locs(iter);
    while flag
        
        if current_index < 2
            flag = 0; start_index(iter) = 1; continue;
        end
        
        if overall_map_smooth(current_index) < threshold_rate && overall_map_smooth(current_index - 1) < threshold_rate
            flag = 0; start_index(iter) = current_index - 1; continue;
        end
        
        current_index = current_index - 1;
    end
    %
    
    % Find end point of the field
    flag = 1;
    current_index = peak_locs(iter);
    while flag
        
        if current_index > map_length - 1
            flag = 0; end_index(iter) = map_length; continue;
        end
        
        if overall_map_smooth(current_index) < threshold_rate && overall_map_smooth(current_index + 1) < threshold_rate
            flag = 0; end_index(iter) = current_index + 1; continue;
        end
        
        current_index = current_index + 1;
    end
    %
    
end



% % Extracted field information

binary_map = zeros(size(overall_map));

for iter = 1 : length(peak_locs)
    binary_map(start_index(iter) : end_index(iter)) = 1;
end

start_index = []; end_index = [];

field_count = 0;
field_index = [];
field_rate = [];
field_size = [];
field_flag = false;

% Find place fields
for iter = 1 : length(binary_map)
    if binary_map(iter) == 1 && field_flag == false
        field_count = field_count + 1;
        field_flag = true;
        start_index(field_count) = iter;
    elseif binary_map(iter) == 0 && field_flag == true
        field_flag = false;
        end_index(field_count) = iter - 1;
    end
end

if field_flag == true
    end_index(field_count) = iter;
end
%


if field_count == 0
    field_count=0; start_index=0; end_index=0; field_size=0; baseline_fr=0; h=0;
    return;
end


% Extract field information
for iter = 1 : field_count
    field_size(iter) = end_index(iter) - start_index(iter) + 1;
    field_rate(iter) = max(overall_map_smooth(start_index(iter) : end_index(iter)));
    field_index(iter) = max(find(overall_map_smooth(start_index(iter) : end_index(iter)) == field_rate(iter))) + start_index(iter) - 1;
    field_center(iter) = (start_index(iter) + end_index(iter)) / 2;
end
%


% Insert nan bins
start_index = start_index + nan_index(1);
end_index = end_index + nan_index(1);
field_index = field_index + nan_index(1);
field_center = field_center + nan_index(1);
%

% %



% % Display

h=[];
% field_boundary_display_2f2;
% field_boundary_function_2f2_epochFiltered;
% %

end
