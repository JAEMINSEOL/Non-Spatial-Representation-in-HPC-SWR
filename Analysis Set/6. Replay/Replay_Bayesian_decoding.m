function [posterior,p_test,R_actual,R_shuffled,v,c] = Replay_Bayesian_decoding(ROOT,Params,thisRip,ReactTable,mapRUNs)
tbinDuration = Params.tbinDuration;
fitting_threshold = Params.fitting_threshold;


RipID = thisRip.ID{1};
thisReact = ReactTable(strcmp(ReactTable.RippleID,RipID),:);
ripple_duration = thisRip.RippleDuration;
nCells = size(unique(thisReact.UnitID),1);
    % time binning

    tbinN = ceil(ripple_duration / tbinDuration);

    %% get N and f(x)
    matRUN = 1; % only CA1

    N = nan(tbinN, nCells);
 mapMat = {}; % tuning curves for all cells (skaggs map 1D)
% Reactivated unit list
if nCells>=3
for UnitNum=1:nCells

    clusterID = thisReact.UnitID{UnitNum};
    findHYPHEN = find(clusterID == '-');
    thisRID = jmnum2str(str2double(clusterID(1, 1:findHYPHEN(1) - 1)),3);
    thisSID = jmnum2str(str2double(clusterID(1, findHYPHEN(1) + 1:findHYPHEN(2) - 1)),2);
    thisTTID = num2str(str2double(clusterID(1, findHYPHEN(2) + 1:findHYPHEN(3) - 1)));
    thisCLID = jmnum2str(str2double(clusterID(1, findHYPHEN(3) + 1:end)),2);

    try
    Spk = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\TT' thisTTID '\parsedSpike_' num2str(str2double(thisCLID)) '.mat']);
%     load([ROOT.Unit '\' clusterID '.mat']);

    % linearized firing map for each reactivated unit
    thisField = table;
    temp=[];
    temp = Spk.t_spk(~Spk.area_spk(:,5) & Spk.correctness_spk(:,1));
    temp = sortrows(temp,1);

    thisField.ts = temp(:,1);
    thisMap = getFieldMaps(clusterID,thisField,'session',ROOT.Raw.Mother,ROOT.Info);

    % get N matrix
    thisSPK = Spk.t_spk(Spk.t_spk>=thisRip.STtime & Spk.t_spk<=thisRip.EDtime);
    thisSPK = thisSPK - thisRip.STtime;
%     thisSPK = thisReact.SpkTime_fromRippleStart(strcmp(clusterID,thisReact.UnitID),:);

    for tbinRUN = 1 : tbinN
        if tbinRUN < tbinN
            thisTimeBin = [0 tbinDuration] + (tbinRUN-1) * tbinDuration;
        else
            thisTimeBin = [0 tbinDuration] + (tbinRUN-1) * tbinDuration;
            thisTimeBin(end) = ripple_duration;
        end
        N(tbinRUN,UnitNum) = length(find(thisSPK >= thisTimeBin(1) & thisSPK < thisTimeBin(2)));
    end

    % get f(x) according to trial information
    if isfield(thisMap, 'skaggsMap1D')
        for i = 1 : size(thisMap.skaggsMap1D,2)
            thisMap.skaggsMap1D{i}(isnan(thisMap.skaggsMap1D{i})) = [];
            mapMat{i}(:,UnitNum) = thisMap.skaggsMap1D{i};
        end
    else
        for mapRUN = 1 : length(pbinN), mapMat{mapRUN}(:,cellRUN)= nan(pbinN(mapRUN),1); end
    end
    end
end

skaggsMap1D = thisMap.skaggsMap1D;
pbinN = [size(skaggsMap1D{1},1),size(skaggsMap1D{2},1) size(skaggsMap1D{3},1) size(skaggsMap1D{4},1) size(skaggsMap1D{5},1)];
if mapRUNs==0
mapRUNs = length(pbinN);
end
occUniform = {};
for mapRUN = 1 : mapRUNs
    occUniform{mapRUN} = repmat([1/pbinN(mapRUN)], [pbinN(mapRUN) 1]);
end

  %% running Bayesian algorithm
                % initialize
                posterior = {};
                numOfLine = []; v = {}; c = {};
                R_actual = []; R_shuffled = {};
                p_test = []; type = {};
                
                tic;
                for mapRUN = 1 : mapRUNs
                    %% get Posterior
                    posterior{mapRUN} = Replay_getP(tbinDuration, occUniform{mapRUN}, N, mapMat{mapRUN});
                    
                    %% line fitting
                    switch Params.lineFittingMethod
                        case 'line finding'
                            % using parfor loop
                            [numOfLine(mapRUN), v{mapRUN}, c{mapRUN}, R_actual(mapRUN)] = Replay_lineFinding(posterior{mapRUN}, fitting_threshold);
                            R_shuffled{mapRUN} = Replay_lineFinding_shuffled(occUniform{mapRUN}, N, mapMat{mapRUN}, tbinDuration, fitting_threshold);
                        case 'linear regression'
                            [v{mapRUN}, R_actual(mapRUN)] = Replay_linearRegression(posterior{mapRUN});
                            R_shuffled{mapRUN} = Replay_linearRegression_shuffled(occUniform{mapRUN}, N, mapMat{mapRUN}, tbinDuration);
                            c{mapRUN} =  v{mapRUN}(1);
                            v{mapRUN} =  v{mapRUN}(2);
                            numOfLine(mapRUN) = 1;
                    end
                    
                    %% test significancy
                    if ~isnan(R_actual(mapRUN))
                        p_test(mapRUN) = length(find(R_shuffled{mapRUN} > R_actual(mapRUN))) / length(R_shuffled{mapRUN});
                        
                        if v{mapRUN}(1) > 0, type{mapRUN} = 'forward';
                        elseif v{mapRUN}(1) < 0, type{mapRUN} = 'reverse';
                        else type{mapRUN} = 'NaN';
                        end
                        
                    elseif isnan(R_actual(mapRUN))
                        p_test(mapRUN) = NaN;
                        
                        v{mapRUN} = NaN;
                        c{mapRUN} = NaN;
                        numOfLine(mapRUN) = NaN;
                        type{mapRUN} = NaN;
                        
                    end                    
                    
                end % for mapRUN
                toc;
else
    disp([thisRip.ID{1} ' has less Place cells'])
end

