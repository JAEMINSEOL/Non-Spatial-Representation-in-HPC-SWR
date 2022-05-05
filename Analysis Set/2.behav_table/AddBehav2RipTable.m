
Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R1'];
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R2'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
RippleList_old = readtable([ROOT.Old '\RipplesTable_CA1.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
Experimenter = {'LSM','JS','SEB'};

RipplesTable = RippleList_old;

fd = dir(ROOT.Old);

RipplesTable_all = table;
thisSID_old = '';

RipplesTable.RippleDuration=(RipplesTable.EDtime-RipplesTable.STtime);
for sid=1:size(RipplesTable,1)

        thisSID = [jmnum2str(RipplesTable.rat(sid),3) '-' jmnum2str(RipplesTable.session(sid),2)];
          load([ROOT.Behav '\' thisSID '.mat']);
          BehavTable = readtable([ROOT.Behav '\' thisSID '.xlsx']);
    %% add to RipplesTable
    Idx = knnsearch(Behav.t,RipplesTable.STtime(sid));
    Idx_t = find(Behav.trial_time(:,1)<RipplesTable.STtime(sid));
    
    
    RipplesTable.PosX(sid) = Behav.x(Idx);
    RipplesTable.PosY(sid) = Behav.y(Idx);
    RipplesTable.PosY_linearized(sid) = Behav.y_linearized(Idx);
    RipplesTable.Vx(sid) = Behav.velocity.Vx(Idx);
    RipplesTable.Vy(sid) = Behav.velocity.Vy(Idx);
    RipplesTable.speed(sid) = Behav.velocity.speed(Idx);
    
    if ~isempty(Idx_t)
        Idx_t = max(Idx_t);
        RipplesTable.trial{sid} = [thisSID '-' jmnum2str(Idx_t,3)];
        RipplesTable.context(sid) = Behav.trial_context(Idx_t);
        RipplesTable.correctness(sid) = Behav.trial_correctness(Idx_t);
        RipplesTable.ambiguity(sid) = Behav.trial_ambiguity(Idx_t);
        RipplesTable.area(sid) = Behav.trial_vector(Idx,2);
        RipplesTable.StartTime_fromTrialEnd(sid) = RipplesTable.STtime(sid) - BehavTable.xEnd(Idx_t);
        RipplesTable.EndTime_fromTrialEnd(sid) = RipplesTable.EDtime(sid) - BehavTable.xEnd(Idx_t);
    else
        RipplesTable.trial{sid} = [thisSID '-000'];
        RipplesTable.StartTime_fromTrialEnd(sid) = nan;
        RipplesTable.EndTime_fromTrialEnd(sid) = nan;
    end


end

RipplesTable = RipplesTable(RipplesTable.StartTime_fromTrialEnd>0,:);

writetable(RipplesTable,[ROOT.Save '\RipplesTable_Behav_' thisRegion '.xlsx']);
