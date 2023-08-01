%% Compare_SFMF

Initial_SWRFilter_common;
warning off

ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip5 = [ROOT.Mother '\Processed Data\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R23_sub'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U2'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];

dir = '';
RegionList = {'CA1','SUB'};

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

RipplesTable.SUB = readtable([ROOT.Save '\RipplesTable_SUB_refCA1_forAnalysis' '.xlsx']);
ReactTable.SUB = readtable([ROOT.Save '\ReactTable_SUB_SUB.xlsx']);
UnitsTable.SUB = readtable([ROOT.Save '\UnitsTable_SUB_forAnalysis_TP.xlsx']);
UnitsTable_field.SUB = readtable([ROOT.Save '\UnitsTable_SUB_field_forAnalysis.xlsx']);

RipplesTable.CA1 = readtable([ROOT.Save '\RipplesTable_CA1_forAnalysis' '.xlsx']);
ReactTable.CA1 = readtable([ROOT.Save '\ReactTable_CA1_CA1.xlsx']);
UnitsTable.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_forAnalysis_TP.xlsx']);
UnitsTable_field.CA1 = readtable([ROOT.Save '\UnitsTable_CA1_field_forAnalysis.xlsx']);

FRMaps.SUB = LoadFRMap(ROOT,UnitsTable.SUB);
FRMaps.SUB_field = LoadFRMap(ROOT,UnitsTable_field.SUB);
FRMaps.CA1 = LoadFRMap(ROOT,UnitsTable.CA1);
FRMaps.CA1_field = LoadFRMap(ROOT,UnitsTable_field.CA1);
%% Fill SI and AvgFR (do once)

%  UnitsTable.CA1(isnan(UnitsTable.CA1.RDI_LScene),:) = [];
%   UnitsTable.SUB(isnan(UnitsTable.SUB.RDI_LScene),:) = [];
% 
% UnitSummary = readtable([ROOT.Info '\ClusterSummary.xlsx']);
% 
% % for u = 1:size(UnitsTable.CA1,1)
% %     idx = find(strcmp(UnitsTable.CA1.ID{u},UnitSummary.ID));
% %     UnitsTable.CA1.onMazeAvgFR(u) = UnitSummary.onMazeAvgFR(idx);
% %     UnitsTable.CA1.SI(u) = UnitSummary.SI(idx);
% % end
% 
% ClusterList = readtable([ROOT.Info '\ClusterList_SWR_SUB.xlsx']);
% for u = 1:size(UnitsTable.SUB,1)
%     idx = find(strcmp(UnitsTable.SUB.ID{u},ClusterList.ID));
%     UnitsTable.SUB.onMazeAvgFR(u) = ClusterList.onMazeAvgFR(idx);
% UnitsTable.SUB.onMazeMaxFR(u) = ClusterList.onMazeMaxFR(idx);
%     UnitsTable.SUB.SI(u) = ClusterList.SI(idx);
% end

%% field COM, Union field
for Reg = 1:numel(RegionList)
    thisRegion = RegionList{Reg};
    
    FRMaps_field = zeros(size(UnitsTable.(thisRegion),1),51);
for uid = 1:size(UnitsTable.(thisRegion),1)
    fidx = find(strncmp(UnitsTable_field.(thisRegion).ID,UnitsTable.(thisRegion).ID{uid},12));
    thisFields = UnitsTable_field.(thisRegion)(fidx,:);

    thisFRMaps = LoadFRMap(ROOT,thisFields);

    for fid = 1:size(thisFields,1)
        thisFieldFR=struct;
        thisFieldFR.numOfSpk1D(1) = 1000;
        thisFieldFR.skaggsMap1D{1} = squeeze(thisFRMaps(1,:,fid));
        [field_count, start_index, end_index, field_size, ~] = field_boundary_function_2f2(thisFieldFR, thisFields.ID{fid});

        FRMaps_field(uid,start_index:end_index) = FRMaps_field(uid,start_index:end_index)+1;

        F = FRMaps.SUB_field(1,:,fidx(fid)); F(isnan(F))=0;
        w = regionprops(true(size(F)), F, 'WeightedCentroid');

        thisFields.COM(fid) = w.WeightedCentroid(1);

        UnitsTable_field.(thisRegion).COM(fidx(fid)) = w.WeightedCentroid(1);
    end

    UnitsTable.(thisRegion).NumField(uid) = size(thisFields,1);
    UnitsTable.(thisRegion).FieldSize_TP(uid) = sum(FRMaps_field(uid,:)>0);
    UnitsTable.(thisRegion).mean_COM(uid) = nanmean(thisFields.COM);

end
end

%% Draw Norm. Population Ratemap

thisFRMaps_Norm=struct;
for Reg = 1:numel(RegionList)
    peaks = [];
    thisRegion = RegionList{Reg};
thisFRMaps_Norm.(thisRegion) = squeeze(FRMaps.(thisRegion)(1,:,:))';
for uid = 1:size(thisFRMaps_Norm.(thisRegion),1)

v = thisFRMaps_Norm.(thisRegion)(uid,:); v(isnan(v))=[];
    x = [1:length(v)];
xq = linspace(1,length(v),51);
vq = interpn(x,v,xq,'linear');

    [m,t] = nanmax(vq);
    peaks(uid,1) = t;
    thisFRMaps_Norm.(thisRegion)(uid,:) = smooth(vq ./ m);
end

% [~,i] = sort(peaks);

[~,i] = sort(UnitsTable.(thisRegion).mean_COM);
thisFRMaps_Norm.(thisRegion) = thisFRMaps_Norm.(thisRegion)(i,:);
end


figure;
% sgtitle('FR peak pos sort')
sgtitle('TP COM pos sort')
subplot(1,2,1)
imagesc(thisFRMaps_Norm.CA1)
title('CA1')

subplot(1,2,2)
imagesc(thisFRMaps_Norm.SUB)
title('SUB')
colormap(jet)

%%
    figure
    x1 = UnitsTable.CA1.FieldSize_TP(UnitsTable.CA1.NumField==1);
    x2 = UnitsTable.CA1.FieldSize_TP(UnitsTable.CA1.NumField>1);
    x3 = UnitsTable.SUB.FieldSize_TP(UnitsTable.SUB.NumField==1);
    x4 = UnitsTable.SUB.FieldSize_TP(UnitsTable.SUB.NumField>1);
x = [x1;x2;x3;x4].*2;
g = [repmat({['CA1-SF (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['CA1-MF (n=' num2str(size(x2,1)) ')']},size(x2,1),1);...
    repmat({['SUB-SF (n=' num2str(size(x3,1)) ')']},size(x3,1),1);repmat({['SUB-MF (n=' num2str(size(x4,1)) ')']},size(x4,1),1)];
boxplot(x,g)

ylabel('field width (cm)')
xlabel('')
ylim([0 100])

%%
    figure
    x1 = UnitsTable.CA1.onMazeAvgFR(UnitsTable.CA1.NumField==1);
    x2 = UnitsTable.CA1.onMazeAvgFR(UnitsTable.CA1.NumField>1);
    x3 = UnitsTable.SUB.onMazeAvgFR(UnitsTable.SUB.NumField==1);
    x4 = UnitsTable.SUB.onMazeAvgFR(UnitsTable.SUB.NumField>1);
x = [x1;x2;x3;x4];
g = [repmat({['CA1-SF (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['CA1-MF (n=' num2str(size(x2,1)) ')']},size(x2,1),1);...
    repmat({['SUB-SF (n=' num2str(size(x3,1)) ')']},size(x3,1),1);repmat({['SUB-MF (n=' num2str(size(x4,1)) ')']},size(x4,1),1)];
boxplot(x,g)

ylabel('Mean firing rate (Hz)')
% xlabel('')
% ylim([0 20])

%%
    figure
    x1 = UnitsTable.CA1.SI(UnitsTable.CA1.NumField==1);
    x2 = UnitsTable.CA1.SI(UnitsTable.CA1.NumField>1);
    x3 = UnitsTable.SUB.SI(UnitsTable.SUB.NumField==1);
    x4 = UnitsTable.SUB.SI(UnitsTable.SUB.NumField>1);
x = [x1;x2;x3;x4];
g = [repmat({['CA1-SF (n=' num2str(size(x1,1)) ')']},size(x1,1),1);repmat({['CA1-MF (n=' num2str(size(x2,1)) ')']},size(x2,1),1);...
    repmat({['SUB-SF (n=' num2str(size(x3,1)) ')']},size(x3,1),1);repmat({['SUB-MF (n=' num2str(size(x4,1)) ')']},size(x4,1),1)];
boxplot(x,g)

ylabel('SI score (bit/spike)')
% xlabel('')
% ylim([0 20])