Initial_SWRFilter_common;

thisRID = 526;
thisRID = jmnum2str(thisRID,3);

for thisSID = 1:3
    thisSID = jmnum2str(thisSID,2);
    for thisTTID = 1:24
        cscID = [thisRID '-' thisSID '-' num2str(thisTTID)];
        % load CSC data
        cscData = loadCSC(cscID, ROOT.Raw.Mother, Params.CSCfileTag, Params.exportMODE, Params.behExtraction,'CSC');
        [eeg_expand,Timestamps_expand] = expandCSC(cscData);
        EEG.Raw = eeg_expand(:);
        
        %%
        % 150-250Hz filtering
        Params.cRange = Params.Ripple;
        
        LocalP = [Params.CSCfileTag '_' num2str(Params.cRange(1)) '-' num2str(Params.cRange(2)) 'filtered'];
        EEG.Filtered = FiltLFP(EEG.Raw,Params.Ripple,Params.Fn,'bandpass');
        % eeg reshape
        colN = ceil(numel(EEG.Filtered)/512);
        eeg_expand_reduced_add = [EEG.Filtered; zeros(colN*512-numel(EEG.Filtered),1)];
        
        
        cscData_filtered = cscData;
        cscData_filtered.eeg = reshape(eeg_expand_reduced_add, 512, colN);
        
        %%
        
        Switch = 'butterSwitch_export';
        filename_tail = [Params.CSCfileTag '_150-250filtered'];
        export_Mat2NlxCSC(cscData_filtered, cscData_filtered, [ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID], filename_tail, Switch)
    end
end