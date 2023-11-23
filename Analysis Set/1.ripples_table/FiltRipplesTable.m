Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Processed '\ripples_mat'];
ROOT.Behav = [ROOT.Processed '\behavior_mat'];
ROOT.React = [ROOT.Processed '\react_mat'];
ROOT.Save = [ROOT.Processed];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';

RipplesTable = readtable([ROOT.Save '\RipplesTable_Behav_' thisRegion '.xlsx']);

crit_time = 20;
%%

RipplesTable_filtered = RipplesTable(RipplesTable.StartTime_fromTrialEnd<crit_time,:);

%%

writetable(RipplesTable_filtered,[ROOT.Save '\RipplesTable_Behav_' thisRegion '_' num2str(crit_time) 's.xlsx'])