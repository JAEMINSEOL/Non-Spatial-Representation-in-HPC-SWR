
Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


Experimenter = {'LSM','SEB','JS'};
thisRegion = 'CA1';


thisSID_p='';

BehavTable_all=table;

for sid = 1:size(SessionList,1)
    if SessionList.include(sid) & ismember(SessionList.experimenter{sid},Experimenter)
        thisSID = [jmnum2str(SessionList.rat(sid),3) '-' jmnum2str(SessionList.session(sid),2)];
          Behav = load([ROOT.Raw.Mother '\rat' thisSID(1:3) '\rat' thisSID '\' 'ParsedPosition.mat']);
        
        

        
        %% Mk BehavTable
        
        BehavTable=table;
        for t=1:Behav.total_trial_number
            BehavTable.TrialID{t} = [thisSID '-' jmnum2str(t,3)];
        end
        
        BehavTable.context = Behav.trial_context;
        BehavTable.choice_side = Behav.trial_side;
        BehavTable.ambiguity = Behav.trial_ambiguity;
        BehavTable.choice_side = Behav.trial_side;
        BehavTable.duration_out = Behav.trial_time(:,5) - Behav.trial_time(:,1);
        BehavTable.duration_in = Behav.trial_time(:,6) - Behav.trial_time(:,5);
        BehavTable.duration_ITI = [Behav.trial_time(2:end,1) - Behav.trial_time(1:end-1,6);0];
        BehavTable.start = Behav.trial_time(:,1);
        BehavTable.s1 = Behav.trial_time(:,2);
        BehavTable.s2 = Behav.trial_time(:,3);
        BehavTable.s3 = Behav.trial_time(:,4);
        BehavTable.reward = Behav.trial_time(:,5);
        BehavTable.end = Behav.trial_time(:,6);
        
        save([ROOT.Behav '\' thisSID '.mat'], 'Behav')
        writetable(BehavTable,[ROOT.Behav '\' thisSID '.xlsx'],'writemode','overwrite');
        
        BehavTable_all =[BehavTable_all;BehavTable];
         save([ROOT.Raw.Mother '\rat' thisSID(1:3) '\rat' thisSID '\' 'ParsedPosition.mat'],'-struct','Behav');
    end
end   

writetable(BehavTable_all,[ROOT.Behav '\BehavTable.xlsx'],'writemode','overwrite');

  
