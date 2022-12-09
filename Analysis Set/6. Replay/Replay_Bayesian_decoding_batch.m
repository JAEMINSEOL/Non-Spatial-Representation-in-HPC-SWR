Initial_SWRFilter_common

ROOT.Unit = [ROOT.Save '\units_mat\U1'];
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '.xlsx']);
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R4'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R4'];
if ~exist(ROOT.Save), mkdir(ROOT.Save); end
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end

rng('shuffle');
gcp = parpool(4);

Params.lineFittingMethod = 'linear regression'; % line finding

%% decoding parameters
Params.tbinDuration = 0.005;
Params.N_threshold = 3;
Params.fitting_threshold = 1.5;
a_test = 0.05;

%% variables for display
max_position = [1.1 1.1 0.1 0.1];
imagePosition = [1600 700 800 300];

ADbituVolts = 0.000000076296274187370727 * 10^6;
Replay_Bayesian_decoding_batch

%% Ripple info
for RipNum = 1:size(RipplesTable,1)
    try
        [posterior,p_test,R_actual,R_shuffled,v,c] = ...
            Replay_Bayesian_decoding(ROOT,Params,RipplesTable,ReactTable,RipNum);
        save([ROOT.Save '\' RipplesTable.ID{RipNum} '.mat'],"c","v",'R_shuffled','R_actual','posterior','p_test')
        RipplesTable.DecodingP_all(RipNum) = p_test(1);
        RipplesTable.DecodingP_Z(RipNum) = p_test(2);
        RipplesTable.DecodingP_P(RipNum) = p_test(3);
        RipplesTable.DecodingP_B(RipNum) = p_test(4);
        RipplesTable.DecodingP_M(RipNum) = p_test(5);
    catch
        RipplesTable.DecodingP_all(RipNum) = nan;
        RipplesTable.DecodingP_Z(RipNum) = nan;
        RipplesTable.DecodingP_P(RipNum) = nan;
        RipplesTable.DecodingP_B(RipNum) = nan;
        RipplesTable.DecodingP_M(RipNum) = nan;
        disp([RipplesTable.ID{RipNum} ' error'])
    end
end

writetable(RipplesTable,[ROOT.Save '\RipplesTable_' thisRegion '_ReplayP.xlsx'])

%% Ripples' replay table
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_ReplayP.xlsx']);
pbinN=ones(1,5);

thisFRMapSCALE=2;

for RipNum = 1:size(RipplesTable,1)
    try
        thisRip = RipplesTable(RipNum,:);
        RipID = RipplesTable.ID{RipNum};
        if exist([ROOT.Save '\' RipID  '.mat'])
            load([ROOT.Save '\' RipID  '.mat'])
            thisRID = jmnum2str(thisRip.rat,3);
            thisSID = jmnum2str(thisRip.session,2);

            Pos = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\ParsedPosition.mat']);
            diverging_point = get_divergingPoint(ROOT.Info, thisRID, thisSID);

            diverging_point = diverging_point*0.23;
            stem_end_index =  (max(Pos.y)-diverging_point)/thisFRMapSCALE;

            %% Display Title
            figure('position',[317,63,923,600],'color','w');
            % title
            subplot(8,6,1)
            title([cell2mat(thisRip.ID) ', ' thisRip.experimenter],'fontsize',15)
            axis off

            % trial info
            subplot(9,6,6)
            if strcmp(thisRip.experimenter,'LSM'), CxtList = {'Zebra','Pebbles','Bamboo','Mountains'};
            elseif strcmp(thisRip.experimenter, 'SEB'), CxtList = {'Dot','Square','Zebra','Pebbles'};
            elseif strcmp(thisRip.experimenter, 'JS'), CxtList = {'Forest','','City',''};
                if thisRip.area==5, thisRip.area=0; end
            end
            cxt = CxtList{thisRip.context};
            if thisRip.correctness, corr='Correct'; else, corr='Wrong'; end
            title(['trial ' (thisRip.trial{1}(end-2:end)) ', ' cxt ', ' corr],'fontsize',10)
            axis off
            sceneName=[{'Overall'},CxtList(:)'];
            %% display reconstructed position probability
            for mapRUN = 1 : length(pbinN)

                temp_position = [0.33 0.4 0.12 0.4] + [0.165 0 0 0] * (mapRUN-2);
                ax = subplot('position', temp_position);

                Replay_display_posterior(posterior{mapRUN}, v{mapRUN}, c{mapRUN}, R_actual(mapRUN), R_shuffled{mapRUN},[], ...
                    p_test(mapRUN), 0,stem_end_index, Params.tbinDuration, sceneName{mapRUN}, ax, a_test,[0 .3 0 .2]);

            end
saveas(gca,[ROOT.Fig '\' cell2mat(thisRip.ID) '.png'])

        close all
        end
    catch
        close all
    end
end