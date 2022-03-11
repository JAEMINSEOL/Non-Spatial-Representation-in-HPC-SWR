
Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);

thisRegion = 'CA1';

RipplesTable_all = load([ROOT.Rip '\R1\RipplesTable.mat']);
RipplesTable = RipplesTable_all.RipplesTable_all;
thisSID_p='';

BehavTable_all=table;

for clRip = 1:size(RipplesTable,1)
    
    RipID = cell2mat(RipplesTable.RippleID(clRip));
    thisSID = RipID(1:6);
    
    %%
    if ~strcmp(thisSID, thisSID_p)
        %% edit Behav structure
        Behav = load([ROOT.Raw.Mother '\rat' thisSID(1:3) '\rat' thisSID '\' 'ParsedPosition.mat']);
        
        Behav.Epoch = csvread([ROOT.Raw.Mother '\rat' thisSID(1:3) '\' ['behaviorEpoch_Rat' thisSID(1:3) '.csv']]);
        Behav.SessionST = Behav.Epoch(str2double(thisSID(5:6)),1);
        Behav.SessionED = Behav.Epoch(str2double(thisSID(5:6)),end);
        
        for i=1:Behav.total_trial_number
            if isempty(find(Behav.trial(:,i))), continue; end
            for j=1:5
                if isempty(Behav.t(find(Behav.trial(:,i) & Behav.area(:,j),1)))
                    Behav.trial_time(i,j) =0;
                else
                    Behav.trial_time(i,j) = Behav.t(find(Behav.trial(:,i) & Behav.area(:,j),1));
                end
            end
            Behav.trial_time(i,6) = Behav.t(find(Behav.trial(:,i),1,'last')+1);
        end
        
        Behav.y_linearized = get_linearized_position(ROOT.Mother,ROOT.Raw.Mother,thisSID);
        Behav.diverging_point = get_divergingPoint(ROOT.Info, thisSID(1:3), thisSID(5:6), 'outbound');
        
        
        Behav.trial_vector(:,1) = interp1(Behav.trial_time(:,1),[1:size(Behav.trial_time,1)]',Behav.t, 'previous','extrap');
        [m,Behav.trial_vector(:,2)] = max(Behav.area,[],2);
        for i=1:size(Behav.trial,1)
            if max(Behav.area(i,:))==0
                Behav.trial_vector(i,2)=0;
            end
        end
        
        %% position interpolation
        
        InTrack = Behav.trial_vector(:,2)~=0;
        Behav.x = interp_pos(Behav.x,InTrack,360);
        Behav.y = interp_pos(Behav.y,InTrack,470);
        
        Behav.y_linearized(~InTrack) = 470;
        [Behav.y_linearized,TF] = fillmissing(Behav.y_linearized,'linear');
        
        %% cal speed
        Behav.velocity=table;
        bin_size = 5;
        v = cal_speed(Behav.t, [Behav.x,Behav.y]*0.02, bin_size);
        
        Behav.velocity.Vx=v(:,1);
        Behav.velocity.Vy=v(:,2);
        Behav.velocity.speed = sqrt(v(:,1).^2+v(:,2).^2);

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
        
        %         writetable(BehavTable,[ROOT.Save '\behavior_mat\' thisSID '.xlsx']);
        save([ROOT.Behav '\' thisSID '.mat'], 'Behav')
        
        BehavTable_all =[BehavTable_all;BehavTable];
    end
    
    
    %% add to RipplesTable
    Idx = knnsearch(Behav.t,RipplesTable.StartTime(clRip));
    Idx_t = find(Behav.trial_time(:,1)<RipplesTable.StartTime(clRip));
    
    
    RipplesTable.PosX(clRip) = Behav.x(Idx);
    RipplesTable.PosY(clRip) = Behav.y(Idx);
    RipplesTable.PosY_linearized(clRip) = Behav.y_linearized(Idx);
    RipplesTable.Vx(clRip) = Behav.velocity.Vx(Idx);
    RipplesTable.Vy(clRip) = Behav.velocity.Vy(Idx);
    RipplesTable.speed(clRip) = Behav.velocity.speed(Idx);
    
    if ~isempty(Idx_t)
        Idx_t = max(Idx_t);
        RipplesTable.trial{clRip} = [thisSID '-' jmnum2str(Idx_t,3)];
        RipplesTable.context(clRip) = Behav.trial_context(Idx_t);
        RipplesTable.correctness(clRip) = Behav.trial_correctness(Idx_t);
        RipplesTable.ambiguity(clRip) = Behav.trial_ambiguity(Idx_t);
        RipplesTable.area(clRip) = Behav.trial_vector(Idx,2);
        RipplesTable.StartTime_fromTrialEnd(clRip) = RipplesTable.StartTime(clRip) - BehavTable.end(Idx_t);
        RipplesTable.EndTime_fromTrialEnd(clRip) = RipplesTable.EndTime(clRip) - BehavTable.end(Idx_t);
    else
        RipplesTable.trial{clRip} = [thisSID '-000'];
        RipplesTable.StartTime_fromTrialEnd(clRip) = nan;
        RipplesTable.EndTime_fromTrialEnd(clRip) = nan;
    end
    
    thisSID_p = thisSID;
    
end

writetable(BehavTable_all,[ROOT.Save '\BehavTable.xlsx']);
writetable(RipplesTable,[ROOT.Save '\RipplesTable_Behav_' thisRegion '.xlsx']);
