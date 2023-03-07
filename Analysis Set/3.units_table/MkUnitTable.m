function UnitsTable = MkUnitTable(ROOT,Behav,Spike,thisCLID,thisRegion,TargetTT,Params)
%%
ClusterList = readtable([ROOT.Info '\ClusterList_SWR_' thisRegion '.xlsx']);
UnitSummary = readtable([ROOT.Info '\ClusterSummary.xlsx']);
UnitsTable=[];


try
    idx = find(strcmp(ClusterList.ID,thisCLID));
    unit=table;
    unit=ClusterList(idx(1),:);
%     unit.UnitID {1}= thisCLID;
%     unit.Region = thisRegion;
%     unit.NumOfSpks = ClusterList.nSpks(idx);
%     unit.SpaInfoScore1D = ClusterList.SI(idx);
%     unit.MeanFR = ClusterList.onMazeAvgFR(idx);
%     unit.PeakFR = ClusterList.onMazeMaxFR(idx);
    
    [c,d,~,~] = CalRDI(thisCLID,ROOT,Behav,Spike);
    unit.ReMap_LScene = c(1);
    unit.ReMap_RScene = c(2);
    unit.ReMap_LR = c(3);

    unit.RDI_LScene = d(1);
    unit.RDI_RScene = d(2);
    unit.RDI_LR = d(3);
    
%     temp = find(strcmp(unit.UnitID, table2cell(UnitSummary(:,1))));
%     if ~isempty(temp)
%         unit.NumOfSpks_all = UnitSummary.x_OfSpks(temp(1));
%         unit = [unit,UnitSummary(temp(1),14:end)];
%         unit.numOfSpk_stem1DOut = str2double(cell2mat(unit.numOfSpk_stem1DOut));
%         unit.onmazeAvgFR_stem1DOut = str2double(cell2mat(unit.onmazeAvgFR_stem1DOut));
%         unit.onmazeMaxFR_stem1DOut = str2double(cell2mat(unit.onmazeMaxFR_stem1DOut));
%         unit.SpaInfoScore_stem1DOut = str2double(cell2mat(unit.SpaInfoScore_stem1DOut));
%     else
%         unit.NumOfSpks_all=nan;
%         unit = [unit,array2table(nan(1,40),'VariableNames',ClusterList(1,14:end).Properties.VariableNames)];
%         unit.cellType = {nan};
%     end
    UnitsTable=[UnitsTable;unit];
    
catch
    disp([unit.ID{1} ' is failed'])
end
end



