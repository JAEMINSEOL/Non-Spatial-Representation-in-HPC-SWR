% y position linearization

% revised from Linearization
% Code by Hyunwoo Lee, 2014-Dec-23
% Modified by Jaemin Seol, 2020-09-04
% Modified by Jaemin Seol for SWR project, 2021-04-23

function y_linearized = get_linearized_position(MotherROOT, RawROOT,clusterID)
%% set ROOTs
SPKmotherROOT = [MotherROOT '\Information Sheet'];
LFPmotherROOT = RawROOT;
%    LFPmotherROOT = [RawROOT '\RawData'];
%% get clusterID info
thisRID = clusterID(1:3);
thisSID = clusterID(5:6);


%% load position data
sessionROOT = [LFPmotherROOT '\rat' thisRID '\rat' thisRID '-' thisSID];
load([sessionROOT '\parsedPosition.mat'],'t','y','x','side','trial','correctness','area','cont','ambiguity');

%% get linearized position
divide_number = 50; % for corner 
bin_size = 1; % for arm part

y_linearized = y;
y_linearized(y==0) = nan;

%% get inpolygon & outbound position data

diverging_point = get_divergingPoint(SPKmotherROOT, thisRID, thisSID, 'outbound');
Boundaries;

pos_outbound = logical([]);
pos_outbound = inpolygon(x, y, xEdge, yEdge) & sum(area(:,1:4),2);

%% Define stem part

pos_stem = logical([]);
pos_stem = inpolygon(x, y, stem_boundary(1,:), stem_boundary(2,:));

%% binning corner area (equal angle)

corner_angle = pi/2 / divide_number;

% Left response side
corner_radius = max(diff(corner_boundary_left(1,:)));
corner_origin = [min(corner_boundary_left(1,:)), max(corner_boundary_left(2,:))];

temp_boundary = [];
temp_boundary(1:2,1) = corner_origin;

for bin_iter = 1 : divide_number
    
    [temp_boundary(1,2), temp_boundary(2,2)] = pol2cart((bin_iter - 1) * corner_angle, corner_radius);
    [temp_boundary(1,3), temp_boundary(2,3)] = pol2cart(bin_iter * corner_angle, corner_radius);
    
    temp_boundary(1,2:3) = corner_origin(1,1) + temp_boundary(1,2:3);
    temp_boundary(2,2:3) = corner_origin(1,2) - temp_boundary(2,2:3);
    
    temp_pos = inpolygon(x, y, temp_boundary(1,:), temp_boundary(2,:));
    
    y_linearized(temp_pos & side(:,1) & pos_outbound) = diverging_point - bin_iter;
end

% Right response side
corner_radius = max(diff(corner_boundary_right(1,:)));
corner_origin = [max(corner_boundary_right(1,:)), max(corner_boundary_right(2,:))];

temp_boundary = [];
temp_boundary(1:2,1) = corner_origin;

for bin_iter = 1 : divide_number
    
    [temp_boundary(1,2), temp_boundary(2,2)] = pol2cart(pi - (bin_iter - 1) * corner_angle, corner_radius);
    [temp_boundary(1,3), temp_boundary(2,3)] = pol2cart(pi - bin_iter * corner_angle, corner_radius);
    
    temp_boundary(1,2:3) = corner_origin(1,1) + temp_boundary(1,2:3);
    temp_boundary(2,2:3) = corner_origin(1,2) - temp_boundary(2,2:3);
    
    temp_pos = inpolygon(x, y, temp_boundary(1,:), temp_boundary(2,:));
    
    y_linearized(temp_pos & side(:,2) & pos_outbound) = diverging_point - bin_iter;
end

corner_end_index = diverging_point - divide_number;


%% binning arm part.

% Left side responses
bin = 0;
for bin_iter = min(corner_boundary_left(1,:)) : bin_size * (-1) : min(xEdge)
    
    bin = bin + 1;
    
    temp_pos = x <= bin_iter & x > bin_iter - bin_size;
    
    y_linearized(temp_pos & side(:,1) & pos_outbound & ~pos_stem) = corner_end_index - bin;
end

% Right side responses
bin = 0;
for bin_iter = max(corner_boundary_right(1,:)) : bin_size : max(xEdge)
    
    bin = bin + 1;
    
    temp_pos = x >= bin_iter & x < bin_iter + bin_size;
    
    y_linearized(temp_pos & side(:,2) & pos_outbound & ~pos_stem) = corner_end_index - bin;
end

arm_end_index = corner_end_index - bin + 1;

end


%

