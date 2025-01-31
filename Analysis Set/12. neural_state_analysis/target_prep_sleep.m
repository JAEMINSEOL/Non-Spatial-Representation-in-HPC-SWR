addpath(genpath('D:\HPC-SWR project\Analysis Program'))


Initial_SWRFilter_common;
warning off
ROOT.Processed = [ROOT.Mother '\Processed Data_sleep'];
ROOT.Save = [ROOT.Processed];

%%
clear RipplesTable

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

sleep_suf_list = {'_pre','_post'};
RegionList = {'CA1','SUB'};
suf = '_sleep';
%% add unit ratio
for reg=1:2

    thisR = RegionList{reg};
    thisRegion0 = thisR;
    thisRegion = thisR;
    thisRegion2 = [thisR '_field'];

    RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_forAnalysis.xlsx']);
    ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '.xlsx']);
    UnitsTable = readtable([ROOT.Save '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
    UnitsTable_field = readtable([ROOT.Save '\UnitsTable_' thisRegion '_field_forAnalysis.xlsx']);
      RipplesTable(RipplesTable.rat==232 & RipplesTable.session==4,:) = [];
      UnitsTable(UnitsTable.rat==232 & UnitsTable.session==4,:) = [];


    UnitsTable.RDI_LR = -UnitsTable.RDI_LR;
    UnitsTable_field.RDI_LR = -UnitsTable_field.RDI_LR;
    for rid=1:size(RipplesTable,1)

        thisReact = ReactTable(strcmp(ReactTable.RippleID,RipplesTable.ID(rid)),:);
        thisUnits=table;
        thisUnits.ID = unique(thisReact.UnitID);
        usz=size(thisUnits,1);

       for uid=1:usz
            thisFields = UnitsTable_field(find(contains(UnitsTable_field.ID,thisUnits.ID{uid})),:);
            if isempty(thisFields), continue; end
            thisUnits.RDI_L_min(uid) = min(thisFields.RDI_LScene);
            thisUnits.RDI_L_max(uid) = max(thisFields.RDI_LScene);
            thisUnits.RDI_R_min(uid) = min(thisFields.RDI_RScene);
            thisUnits.RDI_R_max(uid) = max(thisFields.RDI_RScene);
            thisUnits.RDI_C_min(uid) = min(thisFields.RDI_LR);
            thisUnits.RDI_C_max(uid) = max(thisFields.RDI_LR);

                         idx=find(strcmp(thisUnits.ID{uid},UnitsTable.ID));
             if ~isempty(idx)
                 UnitsTable.RDI_L_min(idx) = thisUnits.RDI_L_min(uid);
                  UnitsTable.RDI_L_max(idx) = thisUnits.RDI_L_max(uid);
                  UnitsTable.RDI_R_min(idx) = thisUnits.RDI_R_min(uid);
                  UnitsTable.RDI_R_max(idx) = thisUnits.RDI_R_max(uid);
                  UnitsTable.RDI_C_min(idx) = thisUnits.RDI_C_min(uid);
                  UnitsTable.RDI_C_max(idx) = thisUnits.RDI_C_max(uid);
             end
                  


       end

        if size(thisUnits,2)<6, continue; end

        p=sum(thisUnits.RDI_L_min<-0.1)/usz; q=sum(thisUnits.RDI_L_max>0.1)/usz;
        if p<q, RipplesTable.Ratio_L(rid) = q; else, RipplesTable.Ratio_L(rid) = 1-p; end   % if >0.5, Bamboo; elseif <0.5, Zebra

        p=sum(thisUnits.RDI_R_min<-0.1)/usz; q=sum(thisUnits.RDI_R_max>0.1)/usz;
        if p<q, RipplesTable.Ratio_R(rid) = q; else, RipplesTable.Ratio_R(rid) = 1-p; end   % if >0.5, Mountain; elseif <0.5, Pebble

        p=sum(thisUnits.RDI_C_min<-0.1)/usz; q=sum(thisUnits.RDI_C_max>0.1)/usz;
        if p<q, RipplesTable.Ratio_C(rid) = q; else, RipplesTable.Ratio_C(rid) = 1-p; end   % if >0.5, Left; elseif <0.5, Right
    end

       temp=[];temp2=[];
temp(:,1) = RipplesTable.nRDI_L_max>=5 & RipplesTable.pBinomDev_L_UV<0.05 & RipplesTable.Ratio_L<0.5;
temp(:,3) = RipplesTable.nRDI_L_max>=5 & RipplesTable.pBinomDev_L_UV<0.05 & RipplesTable.Ratio_L>0.5;
temp(:,2) = RipplesTable.nRDI_R_max>=5 & RipplesTable.pBinomDev_R_UV<0.05 & RipplesTable.Ratio_R<0.5;
temp(:,4) = RipplesTable.nRDI_R_max>=5 & RipplesTable.pBinomDev_R_UV<0.05 & RipplesTable.Ratio_R>0.5;
temp2(:,1) = RipplesTable.nRDI_C_max>=5 & RipplesTable.pBinomDev_C_UV<0.05 & RipplesTable.Ratio_C>0.5;
temp2(:,2) = RipplesTable.nRDI_C_max>=5 & RipplesTable.pBinomDev_C_UV<0.05 & RipplesTable.Ratio_C<0.5;

for c=1:4
    RipplesTable.matchedScn(temp(:,c) & RipplesTable.context==c)=1;
end
RipplesTable.matchedCho((temp2(:,1)) & mod(RipplesTable.context,2)==1)=1;
RipplesTable.matchedCho((temp2(:,2)) & mod(RipplesTable.context,2)==0)=1;


        writetable(RipplesTable,[ROOT.Save '\RipplesTable_' thisRegion0 '_final.xlsx'],'WriteMode','replacefile');
    writetable(ReactTable,[ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '_final.xlsx'],'WriteMode','replacefile');
    writetable(UnitsTable,[ROOT.Save '\UnitsTable_' thisRegion '_final.xlsx'],'WriteMode','replacefile');
    writetable(UnitsTable_field,[ROOT.Save '\UnitsTable_' thisRegion '_field_final.xlsx'],'WriteMode','replacefile');
end

%%
load([ROOT.Mother '\Processed Data\processed_pca.mat'],'Target','Rips','Units')

for reg=1:2

thisR = RegionList{reg};
thisRegion0 = thisR;
thisRegion = thisR;
thisRegion2 = [thisR '_field'];

RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_final.xlsx']);
ReactTable = readtable([ROOT.Save '\ReactTable_' thisRegion '_' thisRegion '_final.xlsx']);
UnitsTable = readtable([ROOT.Save '\UnitsTable_' thisRegion '_final.xlsx']);

for sid = 1:size(SessionList)
    thisRID = SessionList.rat(sid);
    thisSID = SessionList.session(sid);

for slp = 1:2
    thisRip = RipplesTable(RipplesTable.rat==thisRID & RipplesTable.session==thisSID & RipplesTable.sleep_pre_post==slp,:);
    thisUnit = UnitsTable(UnitsTable.rat==thisRID & UnitsTable.session==thisSID,:);
if isempty(thisRip), continue; end
Rips.(['r' jmnum2str(thisRID,3) '_s' jmnum2str(thisSID,2) '_' thisR suf sleep_suf_list{slp}]) = thisRip;
Units.(['r' jmnum2str(thisRID,3) '_s' jmnum2str(thisSID,2) '_' thisR suf sleep_suf_list{slp}]) = thisUnit;
Target.(['r' jmnum2str(thisRID,3) '_s' jmnum2str(thisSID,2) '_' thisR suf sleep_suf_list{slp}])=zeros(size(thisRip,1),size(thisUnit,1));
for rid = 1:size(thisRip,1)
    for uid = 1:size(thisUnit,1)
        tar = find(strcmp(ReactTable.RippleID,thisRip.ID(rid)) & strcmp(ReactTable.UnitID,thisUnit.ID(uid)));
        Target.(['r' jmnum2str(thisRID,3) '_s' jmnum2str(thisSID,2) '_' thisR suf sleep_suf_list{slp}])(rid,uid) = length(tar);
    end
% thisReact = ReactTable(strcmp(ReactTable.RippleID,thisRip.ID(rid)),:);
%              temp= table; thisRUnits=table; thisReact = unique(thisReact.UnitID);
%              for u=1:size(thisReact,1)
%              thisRUnits = [thisRUnits; UnitsTable(strcmp((UnitsTable.ID),thisReact(u)),:)];
%              end
%              Target.(['r' jmnum2str(thisRID,3) '_s' jmnum2str(thisSID,2) '_' thisR suf sleep_suf_list{slp} '_category']) = temp;
end

end
end
end

%%
save([ROOT.Mother '\Processed Data\processed_pca.mat'],'Target','Rips','Units')