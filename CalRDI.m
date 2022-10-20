function [c,d,rand_d,p] = CalRDI(clusterID,ROOT,Behav,Spike)
hyp = find(clusterID=='-');
thisRID = (clusterID(1:hyp(1)-1));
thisSID = (clusterID(hyp(1)+1:hyp(2)-1));
TTID = str2double(clusterID(hyp(2)+1:hyp(3)-1));
UnitID = str2double(clusterID(hyp(3)+1:end));
d=nan(1,3); p=[]; data={};
DoRand = 1;
NumRand = 1000;

diverging_point = get_divergingPoint(ROOT.Info, thisRID, thisSID);
diverging_point = diverging_point*0.23;
Boundaries;
Dv =  (max(Behav.y)-diverging_point)/thisFRMapSCALE;
    
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
            try
            thisFieldMap1D = getFieldMaps(clusterID,thisField,'session',ROOT.Raw.Mother,ROOT.Info);
            [field_count, start_index, end_index, field_size, h] = field_boundary_function_2f2(thisFieldMap1D, clusterID);
            thisFieldMap1D_trial = getFieldMaps(clusterID,thisField,'trial',ROOT.Raw.Mother,ROOT.Info);
            trial_set_all =cat(2,[1:Behav.total_trial_number]',Behav.trial_context,Behav.trial_correctness);
            
              c=nan(1,3);  
            c(1) = nancorrcoef(thisFieldMap1D.skaggsMap1D{2},thisFieldMap1D.skaggsMap1D{3});
            c(2) = nancorrcoef(thisFieldMap1D.skaggsMap1D{4},thisFieldMap1D.skaggsMap1D{5});
            c(3) = nancorrcoef(thisFieldMap1D.skaggsMap_left1D,thisFieldMap1D.skaggsMap_right1D);
            
                
            sample_idx = [];
            for scene_iter = 1 : 4
                sample_idx(:,scene_iter) = trial_set_all(:,2) == scene_iter & trial_set_all(:,3) == 1;
            end
            
            sample_idx = logical(sample_idx);
            sample_idx_cond = {[sample_idx(:,1) sample_idx(:,3)],[sample_idx(:,2) sample_idx(:,4)],[sample_idx(:,1)|sample_idx(:,3) sample_idx(:,2)|sample_idx(:,4)]};
            
            pf=[]; crit_pf = .33;
%                 pf(:,1) = getFRfields(Spk_processed.skaggsMap1D{1,1},crit_pf);
            m=max(45,size(thisFieldMap1D_trial,1));
            pf(:,1) = zeros(1,m);
            if field_count>0
            pf(start_index:end_index)=1;
            end
            catch
                 thisFieldMap1D_trial = getFieldMaps(clusterID,thisField,'trial',ROOT.Raw.Mother,ROOT.Info);
                 m=max(45,size(thisFieldMap1D_trial,1));
                 pf(:,1) = zeros(1,m);
                field_count=0;
                c=nan(1,3);  
            end
            
            if Dv==0, Dv=length(pf); end
            
            rand_d=[];
            for cond_iter = 1 : 3
                if cond_iter==3
                    pf(Dv+1:end,1)=0;
                end
                if sum(pf(:,1))==0, d(1,cond_iter)=nan; rand_d(1:NumRand,cond_iter)=nan; p(1,cond_iter)=nan; continue; end
                id = sample_idx_cond{cond_iter}(:,1);
                sample1 = nanmean(thisFieldMap1D_trial(find(pf(:,1)),id),1);
                id = sample_idx_cond{cond_iter}(:,2);
                sample2 = nanmean(thisFieldMap1D_trial(find(pf(:,1)),id),1);
                
                
                if DoRand
                for rand_iter = 1:NumRand
                    sample_sum = [sample1,sample2];
                    rand_sample = sample_sum(randperm(length(sample_sum)));
                    rand_sample1 = rand_sample(1:length(sample1));
                    rand_sample2 = rand_sample(length(sample1)+1:end);
                    rand_d(rand_iter,cond_iter) = computeCohen_d_noAbs(rand_sample1, rand_sample2);
                end
                end
                
                
                d(1,cond_iter) = computeCohen_d_noAbs(sample1, sample2);
                
%                 [~, p(1,cond_iter)] = ttest(rand_d(:,cond_iter),d(1,cond_iter));
                
                p(1,cond_iter) = min([sum(rand_d(:,cond_iter)<d(1,cond_iter))/NumRand,  sum(rand_d(:,cond_iter)>d(1,cond_iter))/NumRand]);
                
                data{cond_iter*2-1}=sample1;
                data{cond_iter*2}=sample2;
            end
            
        else
           data=NaN; p=NaN;
        end
    else
        data=NaN; p=NaN;
    end
    
    %% 임시
    save([ROOT.Unit '\' clusterID '.mat'],'thisFieldMap1D_trial','d','data','p','rand_d')
    



end
    function c=nancorrcoef(x,y)
x=x(~isnan(x));
y=y(~isnan(y));
if length(x)==length(y)
        c = corrcoef(x,y);
        c=c(2);
else
c=nan;
end
    end
