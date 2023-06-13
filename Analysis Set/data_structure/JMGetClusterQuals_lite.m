function [nSPKS, FRRate,withinREFRACPortion, max_width, max_peak, max_amp, peak_ratio, LRATIO, ISODIST,  LogISIPEAKTIME, valley_slope, valley_proportion,SpaInfoScore,onmazeAvgFR,onmazeMaxFR,offmazeAvgFR] = JMGetClusterQuals_lite(clusterID, ROOT,exper)
sFr = 1 / 32000;
microSEC = 10^6;
mSEC = 10^3;
%it's not correct value. but I don't use the value made using this variables. 2/14/2014
%We used these variables. 2014-Nov.
imROW = 500;
imCOL = 650;
thisFRMapSCALE = 10;
videoSamplingRate = 30;

szDOT = 3;
colTRACE = [.4 .4 .4];
colSPK = [1 0 0];

szFONT = 8;
txtINIX = 1; txtINIY = 9.5; txtADJ = 1;

ISIREFRACTORY = 1; isiWINDOW = 7; isiSCALE = 100; histEDGE = -1:1 / isiSCALE:isiWINDOW;
szLINE = 2;

% define polygon (these values are determined by LSM)
% track type: 1 for short track, 2 for long track, 3 for additional recording (rat561)


close all; %picID = figure('Color', [1 1 1], 'Position', [50 50 800 1000]);

%Parse clusterID
findHYPHEN = find(clusterID == '-');

thisRID = jmnum2str(str2double(clusterID(1, 1:findHYPHEN(1) - 1)),3);
thisSID = jmnum2str(str2double(clusterID(1, findHYPHEN(1) + 1:findHYPHEN(2) - 1)),2);
thisTTID = num2str(str2double(clusterID(1, findHYPHEN(2) + 1:findHYPHEN(3) - 1)));
thisCLID = jmnum2str(str2double(clusterID(1, findHYPHEN(3) + 1:end)),2);


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
if exist(['TT' thisTTID '_beh_SS_' thisCLID '.ntt']) && strcmp(exper,'LSM')
    [thisEpochCLTS, thisEpochCLAP] = Nlx2MatSpike(['TT' thisTTID '_beh_SS_' thisCLID '.ntt'], [1 0 0 0 1], 0, 4, [epochST, epochED]);
        d = dir(['TT' thisTTID '_beh_SS_' thisCLID '.ntt']);
        if str2double(thisCLID)<10
    movefile(d.name, ['TT' thisTTID '_beh_SS_' num2str(str2double(thisCLID)) '.ntt'])
        end
elseif exist(['TT' thisTTID '_beh_SS_' num2str(str2double(thisCLID)) '.ntt'])
    [thisEpochCLTS, thisEpochCLAP] = Nlx2MatSpike(['TT' thisTTID '_beh_SS_' num2str(str2double(thisCLID)) '.ntt'], [1 0 0 0 1], 0, 4, [epochST, epochED]);

end

if strcmp(exper,'JS')
 [~,header] = Nlx2MatSpike(['TT' thisTTID '_beh_SS_' num2str(str2double(thisCLID)) '.ntt'], [1 0 0 0 0], 1, 4, [epochST, epochED]);
 adBitStr = header{contains(header, 'ADBitVolts')};
adBits = sscanf(adBitStr(13:end), '%f');
scaleV = 1 / (adBits(1) * 10^6); %100uV
thisEpochCLAP = thisEpochCLAP ./ scaleV;
else
    thisEpochCLAP = thisEpochCLAP ./ 100;
end
%



nSPKS = size(thisEpochCLTS, 2); %to calculate firing rate


%spike constraint means %1/1/2014 by SB
if nSPKS <= 10
    FRRate = -2;
    onmazeMaxFR = -2;
    onmazeAvgFR = -2;
    offmazeAvgFR = -2;
    offmazeMaxFR = -2;
    withinREFRACPortion = -2;
    max_width = -2;
    max_peak = -2;
    max_amp = -2;
    peak_ratio = -2;
    LRATIO = -2;
    ISODIST = -2;
    SpaInfoScore = -2;
    LogISIPEAKTIME = -2;
    meanVector = [0 0 0 0];
    onmazeMinFR = -2; TMI1 = -2; fft_max_frequency = -2; fft_peak = -2; fft_mean = -2; sparsity = -2; selectivity = -2; fft_selectivity = -2;
    valley_slope = -2; valley_proportion = -2;
    disp([clusterID ' : not sufficient spikes']);
    return;
end

% wholeCLTS = Nlx2MatSpike(['TT',thisTTID,'_beh.ntt'],[1 0 0 0 0], 0, 1, 0);
if exist(['TT',thisTTID,'_beh.ntt'])
wholeCLTS = Nlx2MatSpike(['TT',thisTTID,'_beh.ntt'],[1 0 0 0 0], 0, 1, 0);
elseif exist(['TT',thisTTID,'.ntt'])
    wholeCLTS = Nlx2MatSpike(['TT',thisTTID,'.ntt'],[1 0 0 0 0], 0, 1, 0);
        d = dir(['TT',thisTTID,'.ntt']);
    movefile(d.name, ['TT',thisTTID,'_beh.ntt'])
end
    

ISODIST = nan;
LRATIO = nan;


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
%

% Proportion of valley (area)

start_point = 1;

for iter = peak_point(1) : -1 : 1
    if thisMEANAP(iter, max_channel) < 0
        start_point = iter + 1;
        break;
    end
end


middle_point = peak_point(1);

for iter = valley_point(1) : -1 : 1
    if thisMEANAP(iter, max_channel) > 0
        middle_point = iter;
        break;
    end
end

peak_area = sum(thisMEANAP(start_point : middle_point, max_channel));
valley_area = sum(thisMEANAP(middle_point + 1 : baseline_point - 1, max_channel));

valley_proportion = -valley_area / (peak_area - valley_area);
%


%Inter-spike interval (ISI)
isiHIST = histc(log10(diff(thisEpochCLTS)), histEDGE);
withinREFRACPortion = (sum(diff(thisEpochCLTS) < (ISIREFRACTORY * 10^3)) / length(thisEpochCLTS)) * 100;
LogISIPEAKTIME = (10^histEDGE(min(find(isiHIST == min(max(isiHIST)))))) / 1000;

% load Nlx2MatVT result file
% load('thisPos.mat');

Pos = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\ParsedPosition.mat']);
if str2double(thisCLID) < 10 & exist([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' thisCLID '.mat'])
    d = dir([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' thisCLID '.mat']);
    movefile(d.name, [ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' num2str(str2double(thisCLID)) '.mat'])
end
Spk = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' num2str(str2double(thisCLID)) '.mat']);

if ~isfield(Pos,'cont')
    Pos.cont = Pos.sc;
end

if ~isfield(Pos,'trial_context')
    Pos.cont = Pos.sc;
    Pos.trial_context = max(find(sum(Pos.sc,1)));
end

    [r,c] = find(Pos.cont);
    cxt = zeros(size(Pos.cont,1),1);
    cxt (r) = c;
        [r,c] = find(Spk.cont_spk);
    cxt_spk = zeros(size(Spk.cont_spk,1),1);
    cxt_spk (r) = c;
    


% track type: 1 for short track, 2 for long track, 3 for additional recording (rat561)
diverging_point = get_divergingPoint(ROOT.Info, thisRID, thisSID);
Boundaries;
switch exper
    case 'LSM'

    case 'JS'
        Pos.x1 = Pos.x;
        Pos.x = Pos.y+ rand(length(Pos.y),1);
        Pos.y=Pos.x1;
        Spk.x1 = Spk.x_spk;
        Spk.x_spk = Spk.y_spk + rand(length(Spk.y_spk),1);
        Spk.y_spk = Spk.x1;

    case 'SEB'

end
thisPos = [Pos.t,Pos.x,Pos.y];
thisCLTSforSpatialInfo = [Spk.t_spk,Spk.x_spk,Spk.y_spk];
FRRate = nSPKS / length(thisPos(:,1)) * videoSamplingRate; % mean firing rate

inp = inpolygon(thisPos(:,2), thisPos(:,3), xEdge, yEdge) & logical(cxt);
thisPos = thisPos(inp,:);
ins = inpolygon(thisCLTSforSpatialInfo(:,2), thisCLTSforSpatialInfo(:,3), xEdge, yEdge) & logical(cxt_spk);
thisCLTSforSpatialInfo = thisCLTSforSpatialInfo(ins,:);

nSPKS_in = length(thisCLTSforSpatialInfo(:,1));
nPos_out = Pos.t(~inp);
nSPKs_out = Spk.t_spk(~ins);

thisFRMapSCALE = 2;

[occMap spkMap rawMap skaggsMap] = abmFiringRateMap([thisCLTSforSpatialInfo(:, 1) thisCLTSforSpatialInfo(:, 2) thisCLTSforSpatialInfo(:, 3)], [thisPos(:, 1),[thisPos(:,2:3)]], imROW / thisFRMapSCALE, imCOL / thisFRMapSCALE, thisFRMapSCALE, videoSamplingRate);

SpaInfoScore = GetSpaInfo(occMap, skaggsMap);
onmazeMaxFR = nanmax(nanmax(skaggsMap));
onmazeAvgFR = nanmean(nanmean(skaggsMap));
onmazeMinFR = nanmin(nanmin(skaggsMap));
offmazeAvgFR = length(nSPKs_out) / length(nPos_out)*videoSamplingRate;

fprintf('\n%s is processed\n', clusterID);
clear thisEpochCLTS thisEpochCLAP nvtTS nvtX nvtY nvtHD;
