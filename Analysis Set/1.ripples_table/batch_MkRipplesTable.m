Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
thisRegion = 'CA1';
Exprimenter = ['LSM'];



fd = dir(ROOT.Old);

RipplesTable_all = table;

for i=1:size(SessionList,1)
    
    if SessionList.include(i) && ismember(SessionList.experimenter(i),Exprimenter)
        thisSID = [jmnum2str(SessionList.rat(i),3) '-' jmnum2str(SessionList.session(i),2)];
        load([ROOT.Old '\' thisSID '.mat']);
        
        Rip = MkRipplesTable(ROOT,thisSID,thisRegion);
        
%         RipplesTable_all = [RipplesTable_all; Rip];
        writetable(Rip,[ROOT.Save '\ripples_mat\R1\RipplesTable_' thisSID '_' thisRegion '.xlsx'], 'overwritesheet')
%         save([ROOT.Save '\ripples_mat\R1\' thisSID '_' thisRegion '.mat'], 'Rip')
    end
    
end

%  load([ROOT.Save '\RipplesTable.mat'])


