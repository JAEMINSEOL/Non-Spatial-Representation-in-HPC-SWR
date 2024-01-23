Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Processed '\ripples_mat'];
ROOT.Behav = [ROOT.Processed '\behavior_mat'];
ROOT.React = [ROOT.Processed '\react_mat'];
ROOT.Units = [ROOT.Processed '\units_mat\U2'];
ROOT.Save = [ROOT.Processed];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);


thisRegion = 'SUB';
thisRegion2 = [thisRegion '_field'];
filt_time = 0;

if filt_time==0, suff = ''; else, suff = ['_' num2str(filt_time) 's']; end


ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);
ReactTable = unique([ReactTable(:,1:2)],'rows');

RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis.xlsx']);
% UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion2 '.xlsx']);
UnitsTable_B = readtable([ROOT.Save '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);
UnitsTable_A = readtable([ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis_TP.xlsx']);

UnitsTable = UnitsTable_A;

%%
for uid = 1:size(UnitsTable_A,1)
    thisFields = UnitsTable_B(strncmp(UnitsTable_B.ID, UnitsTable_A.ID(uid),12),:);
    [m1,f1] = max(abs(thisFields.RDI_LScene)); if m1<0.1, f1=nan; end
    [m2,f2] = max(abs(thisFields.RDI_RScene)); if m2<0.1, f2=nan; end
    [m3,f3] = max(abs(thisFields.RDI_LR)); if m3<0.1, f3=nan; end

    if f1==f3 && f2==f3
        UnitsTable_A.FieldHet(uid)=3;
    elseif f1==f3 || f2==f3
        UnitsTable_A.FieldHet(uid)=2;
    elseif f1==f2
        UnitsTable_A.FieldHet(uid)=1;
    elseif ~isnan(mean([f1,f2,f3]))
        UnitsTable_A.FieldHet(uid)=0;
    else
        UnitsTable_A.FieldHet(uid)=nan;
    end
end

 figure
histogram(UnitsTable_A.FieldHet(UnitsTable_A.NumField>1))
%%
writetable(UnitsTable_A,[ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis_TP.xlsx'],'writemode','replacefile')