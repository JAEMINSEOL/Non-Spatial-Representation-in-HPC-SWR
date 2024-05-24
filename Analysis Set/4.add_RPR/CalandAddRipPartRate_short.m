Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Processed '\ripples_mat'];
ROOT.Behav = [ROOT.Processed '\behavior_mat'];
ROOT.React = [ROOT.Processed '\react_mat'];
ROOT.Units = [ROOT.Processed '\units_mat\U2'];
ROOT.Save = [ROOT.Processed];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);


RegionList = {'SUB','CA1'};
for reg=1:2
thisR = RegionList{reg};
% 

thisRegion0 = thisR;
thisRegion = thisR;
thisRegion2 = [thisR '_field'];

% thisRegion = 'SUB';
% thisRegion2 = [thisRegion '_field'];
filt_time = 0;

if filt_time==0, suff = ''; else, suff = ['_' num2str(filt_time) 's']; end


ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);
ReactTable = unique([ReactTable(:,1:2)],'rows');

RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_final.xlsx']);
% UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion2 '.xlsx']);
UnitsTable_B = readtable([ROOT.Save '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);
UnitsTable_A = readtable([ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis_TP.xlsx']);

UnitsTable = UnitsTable_A;
%% Mk RipCountTable
RipCountTable = unique([RipplesTable.rat, RipplesTable.session],'rows');
RipplesTable.pBinomDev_M_UV = min([RipplesTable.pBinomDev_L_UV,RipplesTable.pBinomDev_R_UV,RipplesTable.pBinomDev_C_UV],[],2);
for r=1:size(RipCountTable,1)
    RipCountTable(r,3) = sum(RipplesTable.rat==RipCountTable(r,1) & RipplesTable.session==RipCountTable(r,2));
    RipCountTable(r,4) = sum(RipplesTable.rat==RipCountTable(r,1) & RipplesTable.session==RipCountTable(r,2) & RipplesTable.pBinomDev_M_UV<0.05);
   RipCountTable(r,5) = sum(RipplesTable.rat==RipCountTable(r,1) & RipplesTable.session==RipCountTable(r,2) & RipplesTable.DecodingP_all<0.05);

end




%%
for clUnit = 1:size(UnitsTable,1)
   UnitID= cell2mat(UnitsTable.ID(clUnit));
   UnitsTable.NumField(clUnit) = sum(strncmp(UnitsTable.ID(clUnit),UnitsTable_B.ID,12));
    id = find(cellfun(Params.cellfind(UnitID),(ReactTable.UnitID)));
    thisReactTable = ReactTable(id,:);
    
    for clRip = 1:size(thisReactTable,1)
        RipID = cell2mat(thisReactTable.RippleID(clRip));
        rid = find(cellfun(Params.cellfind(RipID),(RipplesTable.ID)));
        if ~isempty(rid)
            thisReactTable.Ensemble(clRip) = RipplesTable.ensemble(rid);
            thisReactTable.pBinomDev_M_UV(clRip) = RipplesTable.pBinomDev_M_UV(rid);
            thisReactTable.DecodingP_all(clRip) = RipplesTable.DecodingP_all(rid);
        end
    end
    
    if ismember('Ensemble', thisReactTable.Properties.VariableNames)
        thisReactTable(thisReactTable.Ensemble<4,:)=[];
        thisReactTable_NonSp = thisReactTable(thisReactTable.pBinomDev_M_UV<0.05,:);
        thisReactTable_Replay = thisReactTable(thisReactTable.DecodingP_all<0.05,:);

        sid = find(str2double(UnitID(1:3))==RipCountTable(:,1) & str2double(UnitID(5:6))==RipCountTable(:,2));
        
        UnitsTable.RipPartRate_all(clUnit) = size(thisReactTable,1) / RipCountTable(sid,3);
        UnitsTable.RipPartRate_NonSp(clUnit) = size(thisReactTable_NonSp,1) / RipCountTable(sid,4);
        UnitsTable.RipPartRate_N_TR(clUnit) = (size(thisReactTable,1)-size(thisReactTable_NonSp,1)) / ...
            (RipCountTable(sid,3)-RipCountTable(sid,4));

        UnitsTable.RipPartRate_Replay(clUnit) = size(thisReactTable_Replay,1) / RipCountTable(sid,5);

        UnitsTable.RipPartRate_N_R(clUnit) = (size(thisReactTable,1)-size(thisReactTable_Replay,1)) / ...
            (RipCountTable(sid,3)-RipCountTable(sid,5));

    end
    
end
%%
writetable(UnitsTable,[ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis_TP.xlsx'],'writemode','replacefile')
end
%%
figure;
subplot(1,2,1); hold on
cdfplot(UnitsTable.RipPartRate_all(UnitsTable.NumField==1))
cdfplot(UnitsTable.RipPartRate_all(UnitsTable.NumField>1))
title('Ripple Participation Rate (all SWRs)'); legend({'SF','MF'})
subplot(1,2,2); hold on
cdfplot(UnitsTable.RipPartRate_NonSp(UnitsTable.NumField==1))
cdfplot(UnitsTable.RipPartRate_NonSp(UnitsTable.NumField>1))
title('Ripple Participation Rate (task-related reactivations)'); legend({'SF','MF'})
%%
% writetable(UnitsTable,[ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis_TP_240310.xlsx'],'writemode','replacefile')
% writetable(UnitsTable,[ROOT.Save '\UnitsTable_RPR_' thisRegion suff '.xlsx'])
%%

rpr = UnitsTable.RipPartRate_all;

ids=[];
for s = 0:4
    ids(s+1,1) = sum(UnitsTable.Selectivity_LR==s );
    ids(s+1,2) = nanmedian(rpr((UnitsTable.Selectivity_LR==s)));
end
figure
subplot(1,2,1)
pie(ids(:,2))
title('Ripple Part. Rate')
subplot(1,2,2)
pie(ids(:,1))
title('Actual cell ratio')

id = (UnitsTable.Selectivity_LR==4);
[h,p] = ttest(rpr(id)-nanmean(id))


id = (UnitsTable.Selectivity_LR==0);
[h,p] = ttest2(rpr(id),rpr(~id),'Vartype','unequal')

%%
bar([mean(rpr(UnitsTable.MultiVar_C)),mean(rpr(~UnitsTable.MultiVar_C))])
[h,p]=ttest2((rpr(UnitsTable.MultiVar_C)),(rpr(~UnitsTable.MultiVar_C)))

scatter(rpr,max([UnitsTable.RDI_hetero_L,UnitsTable.RDI_hetero_R,UnitsTable.RDI_hetero_C],[],2))
scatter(rpr,UnitsTable.RDI_hetero_SC)
cdfplot(rpr(UnitsTable.MultiVar_SC))

hold on

cdfplot(rpr(~UnitsTable.MultiVar_SC))

% end
