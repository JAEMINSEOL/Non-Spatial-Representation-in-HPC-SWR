
Initial_SWRFilter_common;
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
%%
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
thisRegion = 'CA1';
Experimenter = ['LSM'];



fd = dir(ROOT.Save);

RipplesTable_all = table;

for i=1:size(SessionList,1)
    if SessionList.include(i) && ismember(SessionList.experimenter(i),Experimenter)
        thisSID = [jmnum2str(SessionList.rat(i),3) '-' jmnum2str(SessionList.session(i),2)];
    Recording_region_TT = Recording_region({thisSID},:);
    
    TargetTT = find(cellfun(cellfind(thisRegion),table2array(Recording_region_TT)'));
    
    Ripples = load([ROOT.Save '\' thisSID '.mat']);
    %%
    for t1=1:length(TargetTT)
        t=TargetTT(t1);
        MaxT(t,1) = max(Ripples.Ripples.(['TT' num2str(t)]).index(:,2));
        MaxT(t,2) = max(Ripples.Ripples.(['TT' num2str(t)]).ripples(:,2));
    end
    M = max(MaxT(:,1));
    t = max(MaxT(:,2));
    time_ori = t - M/Params_Ripple.Fs;
    
    RippleArray = zeros(M,size(TargetTT,1));
    
    for t1=1:length(TargetTT)
        t=TargetTT(t1);
        C = Ripples.Ripples.(['TT' num2str(t)]).index;
        for i=1:size(C,1)
            RippleArray(C(i,1):C(i,2),t1) = 1;
        end
    end
    RippleVector = sum(RippleArray,2);
    
    
    n=1;ripples_index=[]; Params.NumTTthreshold=3;
    for i=2:size(RippleVector,1)
        if RippleVector(i-1)<Params.NumTTthreshold && RippleVector(i)>=Params.NumTTthreshold
            ripples_index(n,1)=i;
        elseif RippleVector(i-1)>=Params.NumTTthreshold && RippleVector(i)<Params.NumTTthreshold
            ripples_index(n,2)=i;
            n=n+1;
        end
    end
    ripples_index(ripples_index==0) = 1;
    ripples_index(:,3) = time_ori + ripples_index(:,1)/Params_Ripple.Fs;
    ripples_index(:,4) = time_ori + ripples_index(:,2)/Params_Ripple.Fs;
    for i=1:size(ripples_index,1), ripples_index(i,5) = max(RippleVector(ripples_index(i,1):ripples_index(i,2))); end
    writematrix(ripples_index,[ROOT.Save '\' thisSID '_' thisRegion '.csv'])
    end
end