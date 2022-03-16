Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
ROOT.React = [ROOT.Mother '\Processed Data\react_mat'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';
if filt_time==0, suff = ''; else, suff = ['_' num2str(filt_time) 's']; end

RipplesTable = readtable([ROOT.Save '\RipplesTable_Ensemble_' 'CA1' suff '.xlsx']);
UnitsTable = readtable([ROOT.Save '\UnitsTable_filtered_' thisRegion '.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion suff '.xlsx']);
RipCountTable = readtable([ROOT.Save '\RipplesCountTable_' 'CA1' suff '.xlsx']);
        
%%
for clUnit = 1:size(UnitsTable,1)
   UnitID= cell2mat(UnitsTable.UnitID(clUnit));
    id = find(cellfun(cellfind(UnitID),(ReactTable.UnitID)));
    thisReactTable = ReactTable(id,:);
    
    for clRip = 1:size(thisReactTable,1)
        RipID = cell2mat(thisReactTable.RippleID(clRip));
        rid = find(cellfun(cellfind(RipID),(RipplesTable.RippleID)));
        if ~isempty(rid)
            thisReactTable.Ensemble(clRip) = RipplesTable.Ensemble_all(rid);
        end
    end
    
    if ismember('Ensemble', thisReactTable.Properties.VariableNames)
        thisReactTable(thisReactTable.Ensemble<4,:)=[];
        
        
        sid = find(cellfun(cellfind(UnitID(1:6)),(RipCountTable.SID)));
        
        UnitsTable.RipPartRate_all(clUnit) = sum(thisReactTable.RipCxt~=0) / sum(cell2mat(table2cell(RipCountTable(sid,2:5))));
        UnitsTable.RipPartRate_Z(clUnit) = sum(thisReactTable.RipCxt==1) / RipCountTable.Zebra(sid);
        UnitsTable.RipPartRate_P(clUnit) = sum(thisReactTable.RipCxt==2) / RipCountTable.Pebbles(sid);
        UnitsTable.RipPartRate_B(clUnit) = sum(thisReactTable.RipCxt==3) / RipCountTable.Bamboo(sid);
        UnitsTable.RipPartRate_M(clUnit) = sum(thisReactTable.RipCxt==4) / RipCountTable.Mountains(sid);
        UnitsTable.RipPartRate_L(clUnit) = (sum(thisReactTable.RipCxt==1) + sum(thisReactTable.RipCxt==3)) / (RipCountTable.Zebra(sid) + RipCountTable.Bamboo(sid));
        UnitsTable.RipPartRate_R(clUnit) = (sum(thisReactTable.RipCxt==2) + sum(thisReactTable.RipCxt==4)) / (RipCountTable.Pebbles(sid) + RipCountTable.Mountains(sid));
    end
    
end

writetable(UnitsTable,[ROOT.Save '\UnitsTable_RPR_' thisRegion suff '.xlsx'])
%%
