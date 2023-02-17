function UnitsTable = MkUnitTable_blur(ROOT,Behav,Spike,thisCLID,thisRegion,TargetTT,Params)
%%
ClusterList = readtable([ROOT.Info '\ClusterList_SWR_' thisRegion '_blur.xlsx']);
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
    
    [c,d,~,~] = CalRDI_blur(thisCLID,ROOT,Behav,Spike);
    unit.ReMap_No = c(1);
    unit.ReMap_Lo = c(2);
    unit.ReMap_Hi = c(3);
    unit.ReMap_All = c(4);

    unit.RDI_No = d(1);
    unit.RDI_Lo = d(2);
    unit.RDI_Hi = d(3);
    unit.RDI_All = d(4);
    
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



