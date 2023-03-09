Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data (Blur)'];
ROOT.Rip = [ROOT.Save '\ripples_mat'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];
ROOT.React = [ROOT.Save '\react_mat'];


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';
filt_time = 0;

if filt_time==0, suff = ''; else, suff = ['_' num2str(filt_time) 's']; end

RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion  suff '_forAnalysis.xlsx']);
UnitsTable = readtable([ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion suff '.xlsx']);

%%
for clRip = 1:size(RipplesTable,1)
    RipID = cell2mat(RipplesTable.ID(clRip));
    id = find(cellfun(Params.cellfind(RipID),(ReactTable.RippleID)));
    thisReactTable = ReactTable(id,:);
    
    if ~isempty(thisReactTable)
    for clUnit = 1:size(thisReactTable,1)
        UnitID = cell2mat(thisReactTable.UnitID(clUnit));
        uid = find(cellfun(Params.cellfind(UnitID),(UnitsTable.ID)));
        if ~isempty(uid)
            thisReactTable.RDI_LScene(clUnit) = UnitsTable.RDI_LScene(uid);
            thisReactTable.RDI_RScene(clUnit) = UnitsTable.RDI_RScene(uid);
            thisReactTable.RDI_LR(clUnit) = UnitsTable.RDI_LR(uid);
        else
           thisReactTable.RDI_LScene(clUnit) = nan;
            thisReactTable.RDI_RScene(clUnit) = nan;
            thisReactTable.RDI_LR(clUnit) = nan; 
        end
    end
    RipplesTable.NumPCs(clRip) = size(thisReactTable,1);
    
    RipplesTable.RDI_LScene_Ratio(clRip) = sum(thisReactTable.RDI_LScene>0)/size(thisReactTable,1);
    RipplesTable.RDI_LScene_STD(clRip) = nanstd(thisReactTable.RDI_LScene);
    if isnan(sum(thisReactTable.RDI_LScene)), RipplesTable.RDI_LScene_Ratio(clRip) = nan; end
    
    RipplesTable.RDI_RScene_Ratio(clRip) = sum(thisReactTable.RDI_RScene>0)/size(thisReactTable,1);
    RipplesTable.RDI_RScene_STD(clRip) = nanstd(thisReactTable.RDI_RScene);
    if isnan(sum(thisReactTable.RDI_RScene)), RipplesTable.RDI_RScene_Ratio(clRip) = nan; end
    
    RipplesTable.RDI_LR_Ratio(clRip) = sum(thisReactTable.RDI_LR>0)/size(thisReactTable,1);
    RipplesTable.RDI_LR_STD(clRip) = nanstd(thisReactTable.RDI_LR);
    if isnan(sum(thisReactTable.RDI_LR)), RipplesTable.RDI_LR_Ratio(clRip) = nan; end
    else
        RipplesTable.RDI_LScene_Ratio(clRip) = nan;
        RipplesTable.RDI_RScene_Ratio(clRip) = nan;
        RipplesTable.RDI_LR_Ratio(clRip) = nan;
    
    end
end



cmap = jet(256);
v = rescale(RipplesTable.RDI_LScene_Ratio, 1, 256); % Nifty trick!
numValues = length(RipplesTable.RDI_LScene_Ratio)
markerColors = zeros(numValues, 3)
% Now assign marker colors according to the value of the data.
for k = 1 : numValues
    row = round(v(k));
    markerColors(k, :) = cmap(row, :);
end
% Create the scatter plot.
scatter(x, y, [], markerColors);
grid on;
