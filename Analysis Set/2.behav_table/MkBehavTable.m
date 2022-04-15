
Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
Cluster_List = readtable([ROOT.Info '\ClusterList_SWR_CA1_filtered.xlsx'],'ReadVariableNames',true);

Experimenter = {'LSM','SEB','JS'};
thisRegion = 'CA1';


thisSID_p='';

BehavTable_all=table;

for sid = 1:size(SessionList,1)
    if SessionList.include(sid) & ismember(SessionList.experimenter{sid},Experimenter)
        thisSID = [jmnum2str(SessionList.rat(sid),3) '-' jmnum2str(SessionList.session(sid),2)];
        
        
        
        %% edit Behav structure
        Behav = load([ROOT.Raw.Mother '\rat' thisSID(1:3) '\rat' thisSID '\' 'ParsedPosition.mat']);
        
        if exist([ROOT.Raw.Mother '\rat' thisSID(1:3) '\' ['behaviorEpoch_Rat' thisSID(1:3) '.csv']])
            Behav.Epoch = csvread([ROOT.Raw.Mother '\rat' thisSID(1:3) '\' ['behaviorEpoch_Rat' thisSID(1:3) '.csv']]);
        else
            Behav.Epoch = table2array(readtable([ROOT.Raw.Mother '\rat' thisSID(1:3) '\' ['behaviorEpoch_Rat' thisSID(1:3) '.xlsx']]));
        end
        Behav.SessionST = Behav.Epoch(str2double(thisSID(5:6)),1);
        Behav.SessionED = Behav.Epoch(str2double(thisSID(5:6)),end);
        Behav.total_trial_number = size(Behav.trial,2);
        
        if ~isfield(Behav,'area')
            Behav.area(:,1) = sum(Behav.trial,2);
        Behav.area(:,6) = sum(Behav.trial,2)==0;
        end
        
        for i=1:Behav.total_trial_number
            if isempty(find(Behav.trial(:,i))), continue; end
            for j=1:5
                try
                    if isempty(Behav.t(find(Behav.trial(:,i) & Behav.area(:,j),1)))
                        Behav.trial_time(i,j) =0;
                    else
                        Behav.trial_time(i,j) = Behav.t(find(Behav.trial(:,i) & Behav.area(:,j),1));
                    end
                catch
                     Behav.trial_time(i,1) = Behav.t(find(Behav.trial(:,i),1,'first'));
                end
            end
            Behav.trial_time(i,6) = Behav.t(find(Behav.trial(:,i),1,'last')+1);
        end
        
        try
            Behav.y_linearized = get_linearized_position(ROOT.Mother,ROOT.Raw.Mother,thisSID);
            Behav.diverging_point = get_divergingPoint(ROOT.Info, thisSID(1:3), thisSID(5:6), 'outbound');
        catch
            Behav.y_linearized = Behav.y;
            Behav.diverging_point = floor(min(Behav.y));
        end
        
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
        save([ROOT.Raw.Mother '\rat' thisSID(1:3) '\rat' thisSID '\' 'ParsedPosition.mat'],'-struct','Behav');
    end
end   

writetable(BehavTable_all,[ROOT.Save '\BehavTable.xlsx']);

  