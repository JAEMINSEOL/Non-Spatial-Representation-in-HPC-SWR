Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
ROOT.React = [ROOT.Mother '\Processed Data\react_mat'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'SUB';
thisRegion2 = 'SUB_field';
filt_time = 0;

if filt_time==0, suff = ''; else, suff = ['_' num2str(filt_time) 's']; end


ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);

RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_forAnalysis_RDI.xlsx']);
% UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion2 '.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);

UnitsTable = UnitsTable_A;
%% Mk RipCountTable
RipCountTable = unique([RipplesTable.rat, RipplesTable.session],'rows');

for r=1:size(RipCountTable,1)
    RipCountTable(r,3) = sum(RipplesTable.rat==RipCountTable(r,1) & RipplesTable.session==RipCountTable(r,2));
end

%%
for clUnit = 1:size(UnitsTable,1)
   UnitID= cell2mat(UnitsTable.ID(clUnit));
    id = find(cellfun(Params.cellfind(UnitID),(ReactTable.UnitID)));
    thisReactTable = ReactTable(id,:);
    
    for clRip = 1:size(thisReactTable,1)
        RipID = cell2mat(thisReactTable.RippleID(clRip));
        rid = find(cellfun(Params.cellfind(RipID),(RipplesTable.ID)));
        if ~isempty(rid)
            thisReactTable.Ensemble(clRip) = RipplesTable.ensemble(rid);
        end
    end
    
    if ismember('Ensemble', thisReactTable.Properties.VariableNames)
        thisReactTable(thisReactTable.Ensemble<4,:)=[];
        
        
        sid = find(str2double(UnitID(1:3))==RipCountTable(:,1) & str2double(UnitID(5:6))==RipCountTable(:,2));
        
        UnitsTable.RipPartRate_all(clUnit) = sum(thisReactTable.RipCxt~=0) / RipCountTable(sid,3);
   end
    
end
%%

%%
figure;
rpr = UnitsTable.RipPartRate_all;
scatter(rpr,max([UnitsTable.RDI_hetero_L,UnitsTable.RDI_hetero_R,UnitsTable.RDI_hetero_C],[],2))
scatter(rpr,UnitsTable.RDI_hetero_SC)
cdfplot(rpr(UnitsTable.MultiVar_SC))

hold on

cdfplot(rpr(~UnitsTable.MultiVar_SC))
