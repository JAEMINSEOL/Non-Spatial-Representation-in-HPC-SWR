function [d,data,p] = CalRDI(clusterID,ROOT,Behav,Spike)
hyp = find(clusterID=='-');
RatID = str2double(clusterID(1:hyp(1)-1));
SSID = str2double(clusterID(hyp(1)+1:hyp(2)-1));
TTID = str2double(clusterID(hyp(2)+1:hyp(3)-1));
if TTID==14
    a=1;
end
UnitID = str2double(clusterID(hyp(3)+1:end));
d=nan(1,3); p=[]; data={};

if exist([ROOT.Raw.Map '\rat' clusterID '.mat'])
    Spk_processed = load([ROOT.Raw.Map '\rat' clusterID '.mat']);
    if isfield(Spk_processed,'stem_end_index')
        Dv = Spk_processed.stem_end_index;
    else
        Dv=31;
    end
else
    Dv=0;
end
    
    if isfield(Spike, ['TT' num2str(TTID)])
        if isfield(Spike.(['TT' num2str(TTID)]),['Unit' num2str(UnitID)])
            
            thisSpike = Spike.(['TT' num2str(TTID)]).(['Unit' num2str(UnitID)]);
            
            
            thisField = table;
            thisSpike.y_spk_linearized = interp1(Behav.t,Behav.y_linearized,thisSpike.t_spk);
            
            temp=[];
            
            
            
            temp = thisSpike.t_spk(~thisSpike.area_spk(:,5) & thisSpike.correctness_spk(:,1));
            temp = sortrows(temp,1);
            
            thisField.ts = temp(:,1);
            %         thisField.y_linearized = temp(:,2);
            thisFieldMap1D_trial = getFieldMaps(clusterID,thisField,'trial',ROOT.Raw.Mother,ROOT.Info);
            trial_set_all =cat(2,[1:Behav.total_trial_number]',Behav.trial_context,Behav.trial_correctness);
            
            sample_idx = [];
            for scene_iter = 1 : 4
                sample_idx(:,scene_iter) = trial_set_all(:,2) == scene_iter & trial_set_all(:,3) == 1;
            end
            
            sample_idx = logical(sample_idx);
            sample_idx_cond = {[sample_idx(:,1) sample_idx(:,3)],[sample_idx(:,2) sample_idx(:,4)],[sample_idx(:,1)|sample_idx(:,3) sample_idx(:,2)|sample_idx(:,4)]};
            
            pf=[]; crit_pf = .33;
            %     pf(:,1) = getFRfields(Spk_processed.skaggsMap1D{1,1},crit_pf);
            m=min(45,size(thisFieldMap1D_trial,1));
            pf(:,1) = ones(1,m);
            if Dv==0, Dv=length(pf); end
            
            
            for cond_iter = 1 : 3
                if cond_iter==3
                    pf(Dv+1:end,1)=0;
                end
                id = sample_idx_cond{cond_iter}(:,1);
                sample1 = nanmean(thisFieldMap1D_trial(find(pf(:,1)),id),1);
                id = sample_idx_cond{cond_iter}(:,2);
                sample2 = nanmean(thisFieldMap1D_trial(find(pf(:,1)),id),1);
                
                
                d(1,cond_iter) = computeCohen_d_noAbs(sample1, sample2);
                
                [~, p(1,cond_iter)] = ttest2(sample1,sample2);
                
                data{cond_iter*2-1}=sample1;
                data{cond_iter*2}=sample2;
            end
            
        else
           data=NaN; p=NaN;
        end
    else
        data=NaN; p=NaN;
    end
    
end

