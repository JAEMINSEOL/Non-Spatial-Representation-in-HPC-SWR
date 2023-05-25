Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Mother '\Processed Data\ripples_mat\R4'];
ROOT.Fig = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R13_sub'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U2'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];

dir = '';

if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


%%
thisRegion = 'CA1';
thisRegion2 = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_field_RDIs_UV.xlsx']);
% UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion2 '.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_B;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 1;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.05*Params.Fs;
dur = 0.4*Params.Fs;
thisFRMapSCALE=2;
Params.tbinDuration = 0.005;

filter_ns = 'C';
filter_ns2 = 'LR';

for clRip = 1:size(RipplesTable,1)
   RipID= cell2mat(RipplesTable.ID(clRip));
    id = find(cellfun(Params.cellfind(RipID),(ReactTable.RippleID)));
    thisReactTable = ReactTable(id,:);
    
uniqs = unique(thisReactTable.UnitID); uniqs_id=zeros(size(uniqs,1),1);
for u=1:size(uniqs)
    if ismember(uniqs{u},UnitsTable.ID)
    uniqs_id(u,1)=1;
    end
end
uniqs = uniqs(~uniqs_id,:);
 RipplesTable.nPCs(clRip) = size(uniqs,1);
    
end
RipplesTable_p = RipplesTable;
RipplesTable_p(RipplesTable_p.nPCs<3,:)=[];
writetable(RipplesTable,[ROOT.Save '\RipplesTable_' thisRegion2 '_field_RDIs_UV.xlsx'],'writemode','replacefile');
