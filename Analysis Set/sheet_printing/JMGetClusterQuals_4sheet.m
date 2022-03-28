function [sheet,nSPKS,onmazeAvgFR] = JMGetClusterQuals_4sheet(clusterID, ROOT,sessiontype,exper)

sFr = 1 / 32000;
microSEC = 10^6;
mSEC = 10^3;
%it's not correct value. but I don't use the value made using this variables. 2/14/2014
%We used these variables. 2014-Nov.

videoSamplingRate = 30;

szDOT = 3;
colTRACE = [.4 .4 .4];
colSPK = [1 0 0];

szFONT = 8;
txtINIX = 1; txtINIY = 9.5; txtADJ = 1;

ISIREFRACTORY = 1; isiWINDOW = 7; isiSCALE = 100; histEDGE = -1:1 / isiSCALE:isiWINDOW;
szLINE = 2;


close all; 
% picID = figure('Color', [1 1 1], 'Position', [50 50 800 1000]);

%Parse clusterID
findHYPHEN = find(clusterID == '-');

thisRID = jmnum2str(str2double(clusterID(1, 1:findHYPHEN(1) - 1)),3);
thisSID = jmnum2str(str2double(clusterID(1, findHYPHEN(1) + 1:findHYPHEN(2) - 1)),2);
thisTTID = num2str(str2double(clusterID(1, findHYPHEN(2) + 1:findHYPHEN(3) - 1)));
thisCLID = jmnum2str(str2double(clusterID(1, findHYPHEN(3) + 1:end)),2);

if ~strcmp(exper, 'JS')
    if num2str(thisRID) < 400, track_type = 2;
    elseif num2str(thisRID) > 400 && num2str(thisRID) > 500, track_type = 4;
    elseif num2str(thisRID) > 500, track_type = 3;
    end
    
    if track_type == 1
        xEdge = [330 420 420 490 490 260 260 330 330];
        yEdge = [480 480 270 270 180 180 270 270 480];
    elseif track_type == 2
        xEdge = [330 400 400 480 480 250 250 330 330];
        yEdge = [480 480 160 160 70 70 160 160 480];
    elseif track_type == 4 % for rat415
        xEdge = [325 395 395 475 475 245 245 325 325];
        yEdge = [480 480 190 190 100 100 190 190 480];
    elseif track_type == 3 % for rat561
        xEdge = [325 395 395 475 475 245 245 325 325];
        yEdge = [480 480 200 200 110 110 200 200 490];
    end
else
    
end


%Load Epoch information
cd([ROOT.Raw.Mother '\rat' thisRID]);

if exist(['behaviorEpoch_rat' thisRID '.xlsx'])
    thisEPOCH = xlsread(['behaviorEpoch_rat' thisRID '.xlsx']);
else
    thisEPOCH = csvread(['behaviorEpoch_rat' thisRID '.csv']);
end
    epochST = thisEPOCH(str2double(thisSID),1);
    epochED = thisEPOCH(str2double(thisSID),2);


%Load cluster file
cd([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID]);
if exist(['TT' thisTTID '_beh_SS_' thisCLID '.ntt'])
    [thisEpochCLTS, thisEpochCLAP] = Nlx2MatSpike(['TT' thisTTID '_beh_SS_' thisCLID '.ntt'], [1 0 0 0 1], 0, 4, [epochST, epochED]);
elseif exist(['TT' thisTTID '_beh_SS_' num2str(str2double(thisCLID)) '.ntt'])
    [thisEpochCLTS, thisEpochCLAP] = Nlx2MatSpike(['TT' thisTTID '_beh_SS_' num2str(str2double(thisCLID)) '.ntt'], [1 0 0 0 1], 0, 4, [epochST, epochED]);
end
thisEpochCLAP = thisEpochCLAP ./ 100;

nSPKS = size(thisEpochCLTS, 2); %to calculate firing rate

% random sampling of spikes to save time for waveform display
if nSPKS > 50
   r = randi([1 length(thisEpochCLAP)],100,1);
   thisEpochCLAPrand = thisEpochCLAP(:,:,r);
else
    thisEpochCLAPrand = thisEpochCLAP;
end

% align peak
aligned_wave_mat = [];
for ttRUN = 1:4
    [TT_max_amp, TT_max_ind]= max(mean(squeeze(thisEpochCLAP(:,ttRUN,:)),2));
    for clRUN = 1:size(thisEpochCLTS, 2)  
        new_wave = nan(32,1);
        aligned_wave_mat(:,ttRUN,clRUN) = thisEpochCLAP(:,ttRUN,clRUN);
        [max_amp, max_ind] = max(aligned_wave_mat(:,ttRUN,clRUN));
        
        if max_ind ~= TT_max_ind
            max_diff = max_ind - TT_max_ind;
            
            if max_diff > 0
                new_wave(1:end - max_diff) = aligned_wave_mat(1 + max_diff:end, ttRUN, clRUN);
                aligned_wave_mat(:,ttRUN,clRUN) = new_wave;
                
            else
                new_wave(1 - max_diff:end) = aligned_wave_mat(1:end + max_diff, ttRUN, clRUN);
                aligned_wave_mat(:,ttRUN,clRUN) = new_wave;
            end
        end
    end
end

%Spike width [peak to valley; since spike sometimes doesn't come back to the baseline]
%Spike peaks [peak to valley and from baseline]

MEANmaxAPMat = mean(transpose(squeeze(max(thisEpochCLAP))));
max_channel = min(find(MEANmaxAPMat == max(MEANmaxAPMat)));

thisMEANAP = nanmean(aligned_wave_mat,3);
% Find critical points
peak_point(1) = min(find(thisMEANAP(:, max_channel) == max(thisMEANAP(:, max_channel))));
peak_point(2) = thisMEANAP(peak_point(1), max_channel);

valley_point(1) = min(find(thisMEANAP(:, max_channel) == min(thisMEANAP(peak_point(1) : end, max_channel))));
valley_point(2) = thisMEANAP(valley_point(1), max_channel);

temp = find(thisMEANAP(:, max_channel) >= 0);
if isempty(find(temp >= valley_point(1)))
    baseline_point(1) = 32;
else
    baseline_point(1) = min(temp(find(temp >= valley_point(1))));
end
baseline_point(2) = thisMEANAP(baseline_point(1), max_channel);
%

% Spike peaks, width, slope, peak ratio

max_amp = peak_point(2) - valley_point(2);
max_peak = peak_point(2);
peak_ratio = max_peak / max_amp;

max_width = valley_point(1) - peak_point(1);
max_width = max_width * sFr * microSEC;

valley_slope = max(diff(thisMEANAP(valley_point(1) : baseline_point(1), max_channel)));
valley_slope = valley_slope / max_amp;

ISODIST = nan;
LRATIO = nan;

%Inter-spike interval (ISI)
isiHIST = histc(log10(diff(thisEpochCLTS)), histEDGE);
withinREFRACPortion = (sum(diff(thisEpochCLTS) < (ISIREFRACTORY * 10^3)) / length(thisEpochCLTS)) * 100;
LogISIPEAKTIME = (10^histEDGE(min(find(isiHIST == min(max(isiHIST)))))) / 1000;
%Auto-correlogram
[correlogram corrXlabel] = CrossCorr(transpose(thisEpochCLTS) ./ 100, transpose(thisEpochCLTS) ./ 100, 1, 1000);
correlogram((length(corrXlabel)+1)/2) = 0; % to adjust y limit

Pos = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\ParsedPosition.mat']);
Spk = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' num2str(str2double(thisCLID)) '.mat']);

if ~isfield(Pos,'trial_context')
    Pos.trial_context = max(find(sum(Pos.sc,1)));
end
if ~isfield(Spk,'trial_context')
    Spk.trial_context = Spk.cont_spk;
end

    [r,c] = find(Pos.cont);
    cxt = zeros(size(Pos.cont,1),1);
    cxt (r) = c;
        [r,c] = find(Spk.cont_spk);
    cxt_spk = zeros(size(Spk.cont_spk,1),1);
    cxt_spk (r) = c;
    
switch exper
    case 'LSM'
        if num2str(thisRID) < 400, track_type = 2;
        elseif num2str(thisRID) > 400 && num2str(thisRID) <= 500, track_type = 4;
        elseif num2str(thisRID) > 500, track_type = 3;
        end
        
        if track_type == 1
            xEdge = [330 420 420 490 490 260 260 330 330];
            yEdge = [480 480 270 270 180 180 270 270 480];
        elseif track_type == 2
            xEdge = [330 400 400 480 480 250 250 330 330];
            yEdge = [480 480 160 160 70 70 160 160 480];
        elseif track_type == 4 % for rat415
            xEdge = [325 395 395 475 475 245 245 325 325];
            yEdge = [480 480 190 190 100 100 190 190 480];
        elseif track_type == 3 % for rat561
            xEdge = [325 395 395 475 475 245 245 325 325];
            yEdge = [480 480 200 200 110 110 200 200 490];
        end
        
        imROW = 500;
        imCOL = 650;
    case 'JS'
        xEdge = [0 2 2 0 ];
        yEdge = [1000 1000 0 0];
        Pos.x1 = Pos.x;
        Pos.x = Pos.y+ rand(length(Pos.y),1);
        Pos.y=Pos.x1/10;
        Spk.x1 = Spk.x_spk;
        Spk.x_spk = Spk.y_spk + rand(length(Spk.y_spk),1);
        Spk.y_spk = Spk.x1/10;
        imROW = 800;
        imCOL = 200;
    case 'SEB'
        imROW = 500;
        imCOL = 650;
        xEdge = [300 380 380 430 430 250 250 300 300];
        yEdge = [400 400 200 200 150 150 200 200 400];
        
end
thisPos = [Pos.t,Pos.x,Pos.y];
thisCLTSforSpatialInfo = [Spk.t_spk,Spk.x_spk,Spk.y_spk];

in = inpolygon(thisPos(:,2), thisPos(:,3), xEdge, yEdge) & logical(cxt);
thisPos = thisPos(in,:);
inspk = inpolygon(thisCLTSforSpatialInfo(:,2), thisCLTSforSpatialInfo(:,3), xEdge, yEdge) & logical(cxt_spk);
thisCLTSforSpatialInfo = thisCLTSforSpatialInfo(inspk,:);
FRRate = nSPKS / length(thisPos(:,1)) * videoSamplingRate; % mean firing rate
nSPKS_in = sum(inspk);

thisFRMapSCALE = 10;
[occMap spkMap rawMap skaggsMap] = abmFiringRateMap([thisCLTSforSpatialInfo(:, 1) thisCLTSforSpatialInfo(:, 2) thisCLTSforSpatialInfo(:, 3)], thisPos(:, 1:3), imROW / thisFRMapSCALE, imCOL / thisFRMapSCALE, thisFRMapSCALE, videoSamplingRate);

SpaInfoScore = GetSpaInfo(occMap, skaggsMap);
onmazeMaxFR = nanmax(nanmax(skaggsMap));
onmazeAvgFR = nanmean(nanmean(skaggsMap));
onmazeMinFR = nanmin(nanmin(skaggsMap));

if strcmp(sessiontype,'std4zpbm')
    Title = {'Zebra(L)','Pebbles(R)','Bamboo(L)','Mountain(R)'};
elseif strcmp(sessiontype,'std2zp')
    Title = {'Zebra(L)','Pebbles(R)'};
elseif strcmp(sessiontype,'std2ds')
    Title = {'Dot(L)','Square(R)'};
elseif strcmp(sessiontype,'std4dszp')
    Title = {'Dot(L)','Square(R)','Zebra(L)','Pebbles(R)'};
elseif strcmp(sessiontype,'std2fc')
    Title = {'Forest(L)','City(R)'};
end
UnitPropertySheet_Lv1_exc;

