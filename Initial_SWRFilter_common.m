% Initial_SWRFilter_common

clear all; clc; fclose all;

ROOT.Module = 'D:\Modules';
ROOT.Mother = 'D:\HPC-SWR project';
ROOT.Raw.Mother = 'F:\EPhysRawData\RawData';

ROOT.Raw.Map = [ROOT.Raw.Mother '\map files (outbound) epoch_filtered_smoothed'];
ROOT.Raw.Var = [ROOT.Raw.Mother '\variables for display'];
ROOT.Program = [ROOT.Mother '\Analysis Program'];
ROOT.Info = [ROOT.Mother '\Information Sheet'];
ROOT.Save = [ROOT.Mother '\Processed Data'];
addpath(genpath(ROOT.Program))
addpath(genpath(ROOT.Module))
%% set initial parameters
% cell criterions
Params.crit.frlow=1;
Params.crit.fr=10;
Params.crit.width=300;
Params.crit.si=0.5;

% define frequency
Params.Fs = 2000;  % sampling frequency
Params.Fn = Params.Fs/2;  % Nyquist frequency
Params.F0 = 60/Params.Fn; % notch frequency

% define loadCSC variable
Params.CSCfileTag = 'RateReduced';
Params.exportMODE = 1;
Params.behExtraction = 0;

% switch
Params.noiseFilteringSwitch = 1;
Params.saveSwitch = 1;
Params.exportSwitch = 1;

% define data duration
Params.freqN = 2048;
Params.freqLimit = 450;
Params.freqBin = ceil(Params.freqLimit/Params.Fs * Params.freqN);

Params.noise = [350 450]; %noise Range
Params.Ripple = [150 250]; %SWR range
Params.Theta = [4 10];
Params.low = 20;
Params.cRange = Params.Ripple;

Params_Ripple = SetRippleParams(Params);

cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));
cellfindn = @(string)(@(cell_contents)(strncmp(string,cell_contents,11)));
cellfindn2 = @(string)(@(cell_contents)(strncmp(string,cell_contents,3)));

function Params_Ripple = SetRippleParams(Params)
Params_Ripple.Fs = 2000;                                   % sampling frequency
Params_Ripple.gaussianSTD = 40;                            % gaussianSTD / Fs = moving window size for gaussian smoothing (sec) 0.5ms=1idx
Params_Ripple.thresholdSTD = 4;                            % SD
Params_Ripple.beginthresholdSTD = 2;                           % begin & end SD
Params_Ripple.noiseSTD = 100;                               % 11 SD
Params_Ripple.minDuration = 0.04;                          % sec
Params_Ripple.maxDuration = 0.4;                           % sec
Params_Ripple.groupingInterval = 0.04;                     % sec
Params_Ripple.boundarySTD = 1;                             % SD


if isequal(Params.cRange, Params.noise) %noise range
    Params_Ripple.thresholdSTD = 2;                            % SD
    Params_Ripple.noiseSTD = 11;                               %  not used for noise range
    Params_Ripple.groupingInterval = 2;                     % sec
end
end