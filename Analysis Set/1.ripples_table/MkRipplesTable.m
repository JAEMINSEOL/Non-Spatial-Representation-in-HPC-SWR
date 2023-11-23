function Rip = MkRipplesTable(ROOT,thisSID,thisRegion)

ripples_index = readmatrix([ROOT.Old '\' thisSID '_' thisRegion '.csv']);
Recording_region_TT = Recording_region({thisSID},:);
TargetTT = find(cellfun(cellfind(thisRegion),table2array(Recording_region_TT)'));
EEG = LoadEEGData(ROOT, thisSID, TargetTT,Params,Params_Ripple);

Rip=table;
Rip.Rat(1:size(ripples_index,1)) = str2double(thisSID(1:3));
Rip.Session(1:size(ripples_index,1)) = str2double(thisSID(5:6));
Rip.RipID = [1:size(ripples_index,1)]';
for r=1:size(Rip,1), Rip.Region{r} = thisRegion; end
Rip.RippleStartIndex=ripples_index(:,1);
Rip.RippleEndIndex=ripples_index(:,2);
Rip.RippleStart=ripples_index(:,3);
Rip.RippleEnd=ripples_index(:,4);
Rip.NumTT=ripples_index(:,5);
Rip.RippleGap=([ripples_index(2:end,1);ripples_index(end,2)]-ripples_index(:,2))*(1/Params_Ripple.Fs);

temp = (Rip.RippleGap<Params_Ripple.groupingInterval);
j=0;
for i=1:length(temp)-1
    if temp(i)==1
        if temp(i+1)==1
            j=j+1;
        else
            Rip.RippleEndIndex(i-j)=Rip.RippleEndIndex(i+1);
            Rip.RippleStartIndex(i-j+1:i+1)=0;
            j=0;
        end
    end
end
Rip = removevars(Rip, 'RippleGap');
Rip=Rip(find(Rip.RippleStartIndex),:);

Rip.RippleDuration=( Rip.RippleEndIndex- Rip.RippleStartIndex)*(1/Params_Ripple.Fs);
Rip=Rip(find(Rip.RippleDuration>Params_Ripple.minDuration),:);
Rip=Rip(find(Rip.RippleDuration<=Params_Ripple.maxDuration),:);
Rip.RippleGap = ([Rip.RippleStartIndex(2:end);Rip.RippleEndIndex(end)]-Rip.RippleEndIndex)*(1/Params_Ripple.Fs);

Rip.RipID = [1:size(Rip,1)]';

Rip.NumAllTT(1:size(Rip,1)) = size(TargetTT,1);
Rip = movevars(Rip, 'NumTT', 'After', 'NumAllTT');

for j=1:size(Rip,1)
    thisRippleTable = Rip(j,:);
    Idx = [thisRippleTable.RippleStartIndex:thisRippleTable.RippleEndIndex];
    
    EEG_Prop=struct;
    for t=1:length(TargetTT)
        thisEEG = EEG.(['TT' num2str(TargetTT(t))]);
        [EEG_Prop.Amp(t), EEG_Prop.Freq(t), EEG_Prop.power_ripple(t), EEG_Prop.power_theta(t)] = LFP_Properties(thisEEG,Idx);
    end
    Rip.MaxVoltage(j) = max(EEG_Prop.Amp);
    Rip.MeanVoltage(j) = mean(EEG_Prop.Amp);
    Rip.MaxFreq(j) = max(EEG_Prop.Freq);
    Rip.MeanFreq(j) = mean(EEG_Prop.Freq);
    Rip.RipplePower(j) = mean(EEG_Prop.power_ripple);
    Rip.ThetaPower(j) = mean(EEG_Prop.power_theta);
end