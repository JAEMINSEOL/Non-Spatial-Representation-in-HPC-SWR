function TargetTT = GetTargetTT(ROOT,thisSID,thisRegion,Params,crit)

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'readrownames',true);
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);
% Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_' thisRegion '_filtered.xlsx']);

 Recording_region_TT = Recording_region({thisSID},:);
 thisTT_table = TT_table(TT_table.rat==str2double(thisSID(1:3)) & TT_table.session==str2double(thisSID(5:6)),:);

 TargetTT_p = find(cellfun(Params.cellfindn2(thisRegion),table2array(Recording_region_TT)'));

 if crit
     id = knnsearch(thisTT_table.TT,TargetTT_p);
     TargetTT_p = [TargetTT_p,thisTT_table.NoiseRatio_ITI(id),thisTT_table.NumUnits(id)];
     
     TargetTT = TargetTT_p(TargetTT_p(:,2)<crit,1);

     if isempty(TargetTT), TargetTT = TargetTT_p(find(TargetTT_p(:,2)==min(TargetTT_p(:,2))),:); end
 else
 TargetTT = TargetTT_p;
 end
end
 
 