

for sid=1:2  
    %% add to RipplesTable
    Idx = knnsearch(Behav.t,RipplesTable.StartTime(sid));
    Idx_t = find(Behav.trial_time(:,1)<RipplesTable.StartTime(sid));
    
    
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
        RipplesTable.StartTime_fromTrialEnd(sid) = RipplesTable.StartTime(sid) - BehavTable.end(Idx_t);
        RipplesTable.EndTime_fromTrialEnd(sid) = RipplesTable.EndTime(sid) - BehavTable.end(Idx_t);
    else
        RipplesTable.trial{sid} = [thisSID '-000'];
        RipplesTable.StartTime_fromTrialEnd(sid) = nan;
        RipplesTable.EndTime_fromTrialEnd(sid) = nan;
    end

    
end

RipplesTable = RipplesTable(RipplesTable.StartTime_fromTrialEnd>0,:);

writetable(RipplesTable,[ROOT.Save '\RipplesTable_Behav_' thisRegion '.xlsx']);
