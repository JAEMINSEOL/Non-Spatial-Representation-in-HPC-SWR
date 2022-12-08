Initial_SWRFilter_common

ROOT.Unit = [ROOT.Save '\units_mat\U1'];
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '.xlsx']);
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R4'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R0'];
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
