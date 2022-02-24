
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
        
        %% interpolation
        PosX=Behav.x; PosY=Behav.y; time=Behav.t; PosY_linearized = Behav.y_linearized;
        
        InTrack = Behav.trial_vector(:,2)~=0;
        
        temp = PosX(InTrack);
        x = find(~isnan(temp));
        y = temp(x);
        xq = [1:length(temp)]';
        
        yq = round(interp1(x,y,xq,'linear','extrap'));
        PosX(InTrack) = yq;
        
        temp = PosY(InTrack);
        x = find(~isnan(temp));
        y = temp(x);
        xq = [1:length(temp)]';
        
        yq = round(interp1(x,y,xq,'linear','extrap'));
        PosY(InTrack) = yq;
        
        
        PosY(~InTrack) = 470;
        PosX(~InTrack) = 360;
        PosY_linearized(~InTrack) = 470;
        [PosY_linearized,TF] = fillmissing(PosY_linearized,'linear');
        
        
        Behav.x = PosX; Behav.y = PosY; Behav.y_linearized = PosY_linearized;
        %% get speed
        bin_size = 5;
        bin_range = [1:1:bin_size] - ceil(bin_size/2);
        
        this_speed = nan(length(time),1);
        for t_iter = 1 : length(time)
            temp_range = t_iter + bin_range;
            temp_range(temp_range<1 | temp_range>length(time)) = [];
            if length(temp_range) < ceil(bin_size/2)+1, continue; end % speed not assigned at both ends of position data
            
            x_range = PosX(temp_range(1):temp_range(end));
            x_range(isnan(x_range))=[];
            if isempty(x_range), temp_distance=0;
            else
                temp_distance =  x_range(end) - x_range(1);
            end
            temp_time = time(temp_range(end)) - time(temp_range(1));
            
            this_speed(t_iter,1) = (temp_distance*0.02) / temp_time; % centimeter/sec
        end
        for t_iter = 1 : length(time)
            temp_range = t_iter + bin_range;
            temp_range(temp_range<1 | temp_range>length(time)) = [];
            if length(temp_range) < ceil(bin_size/2)+1, continue; end % speed not assigned at both ends of position data
            
            x_range = PosY(temp_range(1):temp_range(end));
            x_range(isnan(x_range))=[];
            if isempty(x_range), temp_distance=0;
            else
                temp_distance =  x_range(end) - x_range(1);
            end
            temp_time = time(temp_range(end)) - time(temp_range(1));
            
            this_speed(t_iter,2) = (temp_distance*0.2) / temp_time; % centimeter/sec
        end
        
        Behav.velocity=table;
        Behav.velocity.Vx=this_speed(:,1);
        Behav.velocity.Vy=this_speed(:,2);
        Behav.velocity.speed = sqrt(this_speed(:,1).^2+this_speed(:,2).^2);
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
    
    
    %%
    Idx = knnsearch(Behav.t,RipplesTable.StartTime(clRip));
    Idx_t = Behav.trial_vector(Idx,1);
    
    RipplesTable.PosX(clRip) = Behav.x(Idx);
    RipplesTable.PosY(clRip) = Behav.y(Idx);
    RipplesTable.PosY_linearized(clRip) = Behav.y_linearized(Idx);
    RipplesTable.Vx(clRip) = Behav.velocity.Vx(Idx);
    RipplesTable.Vy(clRip) = Behav.velocity.Vy(Idx);
    RipplesTable.speed(clRip) = Behav.velocity.speed(Idx);
    
    if ~isnan(Idx_t)
        RipplesTable.trial{clRip} = [thisSID '-' jmnum2str(Idx_t,3)];
        RipplesTable.context(clRip) = Behav.trial_context(Idx_t);
        RipplesTable.correctness(clRip) = Behav.trial_correctness(Idx_t);
        RipplesTable.ambiguity(clRip) = Behav.trial_ambiguity(Idx_t);
        RipplesTable.area(clRip) = Behav.trial_vector(Idx,2);
    else
        RipplesTable.trial{clRip} = [thisSID '-000'];
    end
    
    thisSID_p = thisSID;
    
end

writetable(BehavTable_all,[ROOT.Save '\BehavTable.xlsx']);
writetable(RipplesTable,[ROOT.Save '\RipplesTable_Behav.xlsx']);
