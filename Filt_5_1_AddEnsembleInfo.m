Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
ROOT.React = [ROOT.Mother '\Processed Data\react_mat'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';

RipplesTable = readtable([ROOT.Save '\RipplesTable_Behav_' thisRegion '.xlsx']);
UnitsTable = readtable([ROOT.Save '\UnitsTable_filtered_' thisRegion '.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '.xlsx']);



for clRip = 1:size(RipplesTable,1)
    RipID= cell2mat(RipplesTable.RippleID(clRip));
    id = find(cellfun(cellfind(RipID),(ReactTable.RippleID)));
    thisReactTable = ReactTable(id,:);
    
    unit=[];
    thisUnitTable =[];
    for clUnit=1:size(thisReactTable,1)
        uid = cell2mat(thisReactTable.UnitID(clUnit));
        unit(clUnit) = str2double(uid(8:9))*10^2+str2double(uid(11:12));
        id = find(cellfun(cellfind(uid),(UnitsTable.UnitID)));
        thisUnitTable = [thisUnitTable; UnitsTable(id(1),:)];
        thisReactTable.SI(clUnit) = UnitsTable.SpaInfoScore1D(id(1));
        thisReactTable.PeakBin(clUnit) = UnitsTable.PeakBin(id(1));
    end
    
    if ~isempty(unit)
        [unq_unit,ia,ic] = unique(unit','rows');
        RipplesTable.Spikes(clRip) = length(unit);
        RipplesTable.Ensemble_all(clRip) = length(unq_unit);
        RipplesTable.Ensemble_PC(clRip) = sum(thisUnitTable.SpaInfoScore1D(ia)>=0.5);
        RipplesTable.Ensemble_OnMazePC(clRip) = sum((thisUnitTable.SpaInfoScore1D(ia)>=0.5) & (thisUnitTable.PeakBin(ia)>1));
        
        ReactTable.Ensemble(id) = RipplesTable.Ensemble_all(clRip);
    end
    
    
end


RipplesTable_Trial = RipplesTable(RipplesTable.context>0,:);

RipplesTable_Ensemble = RipplesTable(RipplesTable.Ensemble_all>=3,:);

writetable(RipplesTable_Ensemble,[ROOT.Save '\RipplesTable_Ensemble_' thisRegion '.xlsx'])
