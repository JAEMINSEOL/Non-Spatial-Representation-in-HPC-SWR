function UnitsTable = MkUnitTable(ROOT,Behav,Spike,thisSID,thisRegion,TargetTT)
%%
UnitSummary = readtable([ROOT.Info '\ClusterSummary_SWR.xlsx']);

fd = dir(ROOT.Raw.Map);
UnitsTable=[];
for i=1:size(fd,1)
    try
        id = fd(i).name;
        if length(id)>5
            if str2double(id(4:6))==str2double(thisSID(1:3)) && str2double(id(8:9))==str2double(thisSID(5:6))
                unit=table;
                k=strfind(id, '-');
                unit.UnitID {1}= [id(4:k(1)-1) '-'  id(k(1)+1:k(2)-1) '-'...
                    jmnum2str(str2double(id(k(2)+1:k(3)-1)),2) '-' id(k(3)+1:k(3)+2)];
                if ismember(str2double(id(k(2)+1:k(3)-1)),TargetTT)
                    
                    CellSpec = load([ROOT.Raw.Map '\' fd(i).name]);
                    
                    if ~isfield(CellSpec,'SpaInfoScore1D'), continue; end
                    
                    unit.Region = thisRegion;
                    unit.NumOfSpks = CellSpec.numOfSpk1D(1);
                    unit.SpaInfoScore1D = CellSpec.SpaInfoScore1D(1);
                    
                    FRbin_temp = CellSpec.skaggsMap1D{1};
                    [unit.PeakFR,unit.PeakBin] = nanmax(FRbin_temp);
                    unit.MeanFR = CellSpec.onmazeAvgFR1D(1);
                    
                    
                    clusterID = id(4:end-4);
                    [d,~,~] = CalRDI(clusterID,ROOT,Behav,Spike);
                    unit.RDI_ZB = d(1);
                    unit.RDI_PM = d(2);
                    unit.RDI_LR = d(3);
                    
                    temp = find(strcmp(unit.UnitID, table2cell(UnitSummary(:,1))));
                    if ~isempty(temp)
                        unit.NumOfSpks_all = UnitSummary.x_OfSpks(temp(1));
                        unit = [unit,UnitSummary(temp(1),14:end)];
                        unit.numOfSpk_stem1DOut = str2double(cell2mat(unit.numOfSpk_stem1DOut));
                        unit.onmazeAvgFR_stem1DOut = str2double(cell2mat(unit.onmazeAvgFR_stem1DOut));
                        unit.onmazeMaxFR_stem1DOut = str2double(cell2mat(unit.onmazeMaxFR_stem1DOut));
                        unit.SpaInfoScore_stem1DOut = str2double(cell2mat(unit.SpaInfoScore_stem1DOut));
                    else
                        unit.NumOfSpks_all=nan;
                        unit = [unit,array2table(nan(1,40),'VariableNames',UnitSummary(1,14:end).Properties.VariableNames)];
                        unit.cellType = {nan};
                    end
                    UnitsTable=[UnitsTable;unit];
                end
            end
        end
    catch
        unit.UnitID{1}
    end
end



