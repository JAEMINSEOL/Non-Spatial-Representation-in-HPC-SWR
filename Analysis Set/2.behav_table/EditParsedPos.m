Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R1'];
ROOT.Unit = [ROOT.Mother '\Processed Data\units_mat\U0'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
ROOT.Save = [ROOT.Mother '\Processed Data'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
Cluster_List = readtable([ROOT.Unit '\ClusterList_SWR_CA1_filtered.xlsx'],'ReadVariableNames',true);

Experimenter = {'JS','LSM','SEB'};
thisRegion = 'CA1';

for sid = 1:size(SessionList,1)
    if SessionList.include(sid) & ismember(SessionList.experimenter{sid},Experimenter)
        thisSID = [jmnum2str(SessionList.rat(sid),3) '-' jmnum2str(SessionList.session(sid),2)];
        
        
        
        %% edit Behav structure
        Behav = load([ROOT.Raw.Mother '\rat' thisSID(1:3) '\rat' thisSID '\' 'ParsedPosition.mat']);
        
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
        
        if strcmp(SessionList.experimenter{sid},'SEB')
            t=1; cx=360; cy =470;
            while sum(Behav.trial(t,:))==0
                t=t+1;
            end
            for i=t:size(Behav.trial,1)
                if sum(Behav.trial(i,:))==0
                    if Behav.y(i)>370
                        Behav.trial(i,:) = 0;
                        Behav.area(i,1:5)=0;
                    else
                        Behav.trial(i,:) = Behav.trial(i-1,:);
                        Behav.area(i,5)=1;
                    end
                else
                    if Behav.y(i)<=220
                        Behav.area(i,4)=1;
                    else
                        Behav.area(i,1)=1;
                    end
                end
            end
        elseif strcmp(SessionList.experimenter{sid},'JS')
            cx=1000; cy =1;
            elseif strcmp(SessionList.experimenter{sid},'LSM')
            cx=360; cy =470;
        end
        if ~isfield(Behav,'area')
            Behav.area(:,1) = sum(Behav.trial,2);
            Behav.area(:,5) = sum(Behav.trial,2)==0;
        elseif isempty(Behav.area)
            Behav.area(:,1) = sum(Behav.trial,2);
            Behav.area(:,5) = sum(Behav.trial,2)==0;
        end
        if ~isfield(Behav,'trial_context')
            Behav.trial_context=[];
            for t=1:size(Behav.trial,2)
                i= min(find(Behav.trial(:,t)));
                if i>0
                    Behav.trial_context(t,1) = find(Behav.sc(i,:));
                else
                    Behav.trial_context(t,1) = 0;
                end
            end
        end
        if ~isfield(Behav,'trial_side')
            Behav.trial_side=[];
            for t=1:size(Behav.trial,2)
                i= min(find(Behav.trial(:,t)));
                if i>0
                    Behav.trial_side(t,1) = find(Behav.side(i,:));
                else
                    Behav.trial_side(t,1) = 0;
                end
            end
        end
        if ~isfield(Behav,'trial_ambiguity')
            Behav.trial_ambiguity=[];
            for t=1:size(Behav.trial,2)
                i= min(find(Behav.trial(:,t)));
                if i>0
                    Behav.trial_ambiguity(t,1) = find(Behav.amb(i,:));
                else
                    Behav.trial_ambiguity(t,1) = 0;
                end
            end
        end
        
        
        
        for i=1:Behav.total_trial_number
            if isempty(find(Behav.trial(:,i))), Behav.trial_time(i,:) =0; continue; end
            for j=1:5
                if isempty(Behav.t(find(Behav.trial(:,i) & Behav.area(:,j),1)))
                    Behav.trial_time(i,j) =0;
                else
                    Behav.trial_time(i,j) = Behav.t(find(Behav.trial(:,i) & Behav.area(:,j),1));
                end
                
            end
            Behav.trial_time(i,6) = Behav.t(find(Behav.trial(:,i),1,'last')+1);
        end
        Behav.trial_time(Behav.trial_time(:,1)==0,:) = nan;
        for i=1:6
            Behav.trial_time(:,i) = fillmissing(Behav.trial_time(:,i),'linear','SamplePoints',[1:size(Behav.trial_time,1)]);
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
        
        Behav.x = interp_pos(Behav.x,InTrack,cx);
        Behav.y = interp_pos(Behav.y,InTrack,cy);
        
        Behav.y_linearized(~InTrack) = cy;
        [Behav.y_linearized,TF] = fillmissing(Behav.y_linearized,'linear');
        
        %% cal speed
        Behav.velocity=table;
        bin_size = 5;
        v = cal_speed(Behav.t, [Behav.x,Behav.y]*0.02, bin_size);
        
        Behav.velocity.Vx=v(:,1);
        Behav.velocity.Vy=v(:,2);
        Behav.velocity.speed = sqrt(v(:,1).^2+v(:,2).^2);
        
        save([ROOT.Raw.Mother '\rat' thisSID(1:3) '\rat' thisSID '\' 'ParsedPosition.mat'],'-struct','Behav');
        disp([thisSID ' is finished!'])
    end
end
