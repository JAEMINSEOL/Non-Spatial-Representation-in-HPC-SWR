Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
ROOT.React = [ROOT.Mother '\Processed Data\react_mat'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';
filt_time = 0;

if filt_time==0, suff = ''; else, suff = ['_' num2str(filt_time) 's']; end

RipplesTable = readtable([ROOT.Save '\ripples_mat\R3\RipplesTable_' thisRegion suff '_forAnalysis.xlsx']);
UnitsTable = readtable([ROOT.Save '\units_mat\U1\UnitsTable_CA1_forAnalysis.xlsx']);

thisSID_p='';

BehavTable_all=table;
        ReactTable=table;
        

%%
for clRip = 1:size(RipplesTable,1)
    
    RipID = cell2mat(RipplesTable.ID(clRip));
    thisSID = RipID(1:6);
    thisRipple = RipplesTable(clRip,:);
    
    %%
    if ~strcmp(thisSID, thisSID_p)
        Recording_region_TT = Recording_region({thisSID},:);
        if strcmp(thisRegion, 'CA3')
            TargetTT = find(cellfun(Params.cellfindn2(thisRegion),table2array(Recording_region_TT)'));
        else
            TargetTT = find(cellfun(Params.cellfind(thisRegion),table2array(Recording_region_TT)'));
        end
        
        Spike = LoadSpikeData(ROOT, thisSID, [1:24],Params.cellfindn);
    end
    %%
    for clUnit = 1:size(UnitsTable,1)
        clusterID = cell2mat(UnitsTable.ID(clUnit));
        if strcmp(clusterID(1:6),thisSID)
            thisSpike = Spike.(['TT' num2str(str2double(clusterID(8:9)))]).(['Unit' num2str(str2double(clusterID(11:12)))]);
            SpkTime = thisSpike.t_spk;
            id = SpkTime>=thisRipple.STtime & SpkTime<=thisRipple.EDtime;
            
            if ~isempty(SpkTime(id))
                thisSpkTime = SpkTime(id);
                for spk = 1:size(thisSpkTime,1)
                    react = table;
                    react.RippleID = RipID;
                    react.UnitID = clusterID;
                    react.SpkTime = thisSpkTime(spk);
                    react.SpkTime_fromRippleStart = thisSpkTime(spk) - thisRipple.STtime;
                    react.RipCxt = thisRipple.context;
%                     react.RDI_LScene = UnitsTable.RDI_ZB(clUnit);
%                     react.RDI_PM = UnitsTable.RDI_PM(clUnit);
%                     react.RDI_LR = UnitsTable.RDI_LR(clUnit);
                end
                ReactTable = [ReactTable; react];
            end
        end
    end
    
    thisSID_p = thisSID;
end

writetable(ReactTable,[ROOT.Save '\ReactTable_' thisRegion suff '.xlsx']);
