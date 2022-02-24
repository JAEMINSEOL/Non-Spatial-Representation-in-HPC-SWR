function RateReducer(thisRID, thisSID, RawROOT, SaveROOT)
for thisCSCID = 1:24
    
cscID = [thisRID '-' thisSID '-' num2str(thisCSCID)];
Params.CSCfileTag = 'RateReduced';
Params.exportMODE = 1;
Params.behExtraction = 0;

cscData_raw = loadCSC(cscID, RawROOT, '', Params.exportMODE, Params.behExtraction);
% cscData_reduced = loadCSC(cscID, RawROOT, 'RateReduced', Params.exportMODE, Params.behExtraction);
% cscData_reduced_Copy = loadCSC(cscID, 'D:\HPC-SWR project\temp', 'RateReduced', Params.exportMODE, Params.behExtraction);
if cscData_raw.Timestamps==0
    disp([cscID ' is empty!']);
    continue;
end

cscData_reduced = cscData_raw;

% eeg, timestamp reducing
[eeg_expand_raw,Timestamps_expand_raw] = expandCSC(cscData_raw);
eeg_expand_reduced = eeg_expand_raw(1:16:end);
Timestamps_expand_reduced = Timestamps_expand_raw(1:16:end);

% eeg reshape
colN = ceil(numel(eeg_expand_reduced)/512);
eeg_expand_reduced_add = [eeg_expand_reduced; zeros(colN*512-numel(eeg_expand_reduced),1)];
cscData_reduced.eeg = reshape(eeg_expand_reduced_add, 512, colN);

% reducing other parameters
cscData_reduced.Timestamps = cscData_raw.Timestamps(1:16:end);
cscData_reduced.SampleFrequencies = floor(cscData_raw.SampleFrequencies(1:16:end) / 16);
cscData_reduced.ChannelNumbers = cscData_raw.ChannelNumbers(1:16:end);
cscData_reduced.NumberOfValidSamples = cscData_raw.NumberOfValidSamples(1:16:end);

%header
cscData_reduced.Header(2) = {['## FileName: ' [SaveROOT '\rat' thisRID '\rat' thisRID '-' thisSID] '\' ['CSC' cscData_raw.thisCSCID '_RateReduced.ncs']]};
c=clock;
cscData_reduced.Header(3) = {['## Time Opened: (m/d/y): ' num2str(c(2)) '/' num2str(c(3)) '/' num2str(c(1))...
    ' At Time: ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6))]};
cscData_reduced.Header(4) = cscData_raw.Header(5);
cscData_reduced.Header(5) = {['-NLX_Base_Class_Name']};
cscData_reduced.Header(6) = {['-NLX_Base_Class_Type CscAcqEnt']};
cscData_reduced.Header(7) = cscData_raw.Header(9);
cscData_reduced.Header(8) = cscData_raw.Header(18);
cscData_reduced.Header(9) = {'-ADGain'};
cscData_reduced.Header(10) = {'-AmpGain'};
cscData_reduced.Header(11) = {'-AmpLowCut'};
cscData_reduced.Header(12) = {'-AmpHiCut'};
cscData_reduced.Header(13) = {'-SubSamplingInterleave'};
cscData_reduced.Header(14) = {['-SamplingFrequency ' num2str(cscData_raw.SampleFrequencies(1) / 16)]};
cscData_reduced.Header(15) = cscData_raw.Header(15);
cscData_reduced.Header(16) = cscData_raw.Header(14);
cscData_reduced.Header(17:end) = [];

% export
Switch = 'butterSwitch_export';
filename_tail = 'RateReduced';

cd(SaveROOT)
if ~exist(['rat' thisRID]), mkdir(['rat' thisRID]); end
cd(['rat' thisRID])
if ~exist(['rat' thisRID '-' thisSID]), mkdir(['rat' thisRID '-' thisSID]); end

export_Mat2NlxCSC(cscData_reduced, cscData_reduced, [SaveROOT '\rat' thisRID '\rat' thisRID '-' thisSID], filename_tail, Switch)
end
