Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
ROOT.React = [ROOT.Mother '\Processed Data\react_mat'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';

RipplesTable = readtable([ROOT.Save '\RipplesTable_Behav_' thisRegion '.xlsx']);

crit_time = 20;
%%

RipplesTable_filtered = RipplesTable(RipplesTable.StartTime_fromTrialEnd<crit_time,:);

%%

writetable(RipplesTable_filtered,[ROOT.Save '\RipplesTable_Behav_' thisRegion '_' num2str(crit_time) 's.xlsx'])