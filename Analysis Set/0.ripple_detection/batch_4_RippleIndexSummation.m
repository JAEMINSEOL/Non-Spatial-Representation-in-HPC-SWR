
Initial_SWRFilter_common;
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
%%
Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);
thisRegion = 'CA1';
Experimenter = {'LSM','JS','SEB'};



fd = dir(ROOT.Save);

RipplesTable_all = table;

for sid=1:size(SessionList,1)
    if SessionList.include(sid) & ismember(SessionList.experimenter{sid},Experimenter) 
        thisSID = [jmnum2str(SessionList.rat(sid),3) '-' jmnum2str(SessionList.session(sid),2)];
    Recording_region_TT = Recording_region({thisSID},:);
    
    TargetTT = find(cellfun(cellfind(thisRegion),table2array(Recording_region_TT)'));
    if size(TargetTT,1)<2, continue; end
    load([ROOT.Save '\' thisSID '.mat']);
    %%
    for t1=1:length(TargetTT)
        try
        t=TargetTT(t1);
        MaxT(t,1) = max(Ripples.(['TT' num2str(t)]).index(:,2));
        MaxT(t,2) = max(Ripples.(['TT' num2str(t)]).ripples(:,2));
        end
    end
    M = max(MaxT(:,1));
    t = max(MaxT(:,2));
    time_ori = t - M/Params_Ripple.Fs;
    
    RippleArray = zeros(M,size(TargetTT,1));
    
    for t1=1:length(TargetTT)
        try
        t=TargetTT(t1);
        C = Ripples.(['TT' num2str(t)]).index;
        for i=1:size(C,1)
            RippleArray(C(i,1):C(i,2),t1) = 1;
        end
        end
    end
    RippleVector = sum(RippleArray,2);
    
    
    n=1;ripples_index=[]; Params.NumTTthreshold=min(3,length(TargetTT));
    for i=2:size(RippleVector,1)
        if RippleVector(i-1)<Params.NumTTthreshold && RippleVector(i)>=Params.NumTTthreshold
            ripples_index(n,1)=i;
        elseif RippleVector(i-1)>=Params.NumTTthreshold && RippleVector(i)<Params.NumTTthreshold
            ripples_index(n,2)=i;
            n=n+1;
        end
    end
    
    % grouping
    ripples2=[]; ripples_index2=[];
    i=1;
    while i<size(ripples_index,1)
        j=i+1;
        while ripples_index(j,1) - ripples_index(j-1,2) <= Params_Ripple.groupingInterval*Params_Ripple.Fs
            j=j+1;
            if j>size(ripples_index,1), continue; end
        end
        ripples_index2 = [ripples_index2;[ripples_index(i,1),ripples_index(j-1,2)]];
        i=j;
    end
    if ripples_index(end,2)~=ripples_index2(end,2)
        ripples_index2 = [ripples_index2 ; ripples_index(end,1:2)];
    end
    ripples_index= ripples_index2;
    
    ripples_index(ripples_index(:,2)-ripples_index(:,1) < Params_Ripple.minDuration*Params_Ripple.Fs, :) = [];
    if ripples_index(end,2)<ripples_index(end,1), ripples_index(end,2) = size(RippleVector,1); end
    ripples_index(ripples_index==0) = 1;
    
    %%
    ripples_index(:,3) = time_ori + ripples_index(:,1)/Params_Ripple.Fs;
    ripples_index(:,4) = time_ori + ripples_index(:,2)/Params_Ripple.Fs;
    for i=1:size(ripples_index,1), ripples_index(i,5) = max(RippleVector(ripples_index(i,1):ripples_index(i,2))); end
    
    writematrix(ripples_index,[ROOT.Save '\' thisSID '_' thisRegion '.csv'],'WriteMode', 'overwrite')
    disp([thisSID '_' thisRegion '.csv is saved!'])
    end
end
