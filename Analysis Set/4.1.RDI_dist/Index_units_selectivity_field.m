% Index_units_selectivity

Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Units2 = [ROOT.Mother '\Processed Data\units_mat\U2'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
if ~exist(ROOT.Units2), mkdir(ROOT.Units2); end
dir = '';
thisRSID_old = '';


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


%%
thisRegion = 'CA1';
thisRegion2 = 'CA1_field';
% RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_forAnalysis_RDI.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);
thisFRMapSCALE=2;

Experimenter = {'LSM'};

for uid = 1:size(UnitsTable,1)
    thisUnit = UnitsTable(uid,:);

    UnitID = thisUnit.ID{1};
    [thisRID,thisSID,thisTTID,thisCLID, thisFLID] = parsing_clusterID(UnitID,2);
    [thisRIDn,thisSIDn,thisTTIDn,thisCLIDn, thisFLIDn] = parsing_clusterID(UnitID,8);
    thisRSID =[thisRID '-' thisSID];
    if ~ismember(thisUnit.experimenter,Experimenter), continue; end

    if ~strcmp(thisRSID,thisRSID_old)
        Pos = load([ROOT.Raw.Mother '\rat' thisRID '\rat' thisRID '-' thisSID '\ParsedPosition.mat']);
        diverging_point = get_divergingPoint(ROOT.Info, thisRID, thisSID);

        diverging_point = diverging_point*0.23;
        stem_end_index =  (max(Pos.y)-diverging_point)/thisFRMapSCALE;

        Recording_region_TT = Recording_region({thisRSID },:);

        TargetTT = find(cellfun(Params.cellfindn2(thisRegion),table2array(Recording_region_TT)'));

        if strcmp(thisRegion2,'SUB_field')
            Spike = LoadSpikeData_field(ROOT, thisRSID, TargetTT,Params.cellfindn);
        else
            Spike = LoadSpikeData(ROOT, thisRSID, TargetTT,Params.cellfindn);
        end
        load([ROOT.Behav '\' thisRSID '.mat'])

        thisRSID_old = thisRSID;
    end

    fids = find(UnitsTable_B.rat==thisRIDn & UnitsTable_B.session==thisSIDn & UnitsTable_B.TT==thisTTIDn & UnitsTable_B.AvgPeaktoValley==thisUnit.AvgPeaktoValley);
    %%
    thisUnit_f = UnitsTable_B(fids,:);
    %
    %     UnitsTable.Selectivity_LScene(uid) = Cell_Cat(thisUnit_f,stem_end_index,'LScene');
    %     UnitsTable.Selectivity_RScene(uid) = Cell_Cat(thisUnit_f,stem_end_index,'RScene');
    %     UnitsTable.Selectivity_LR(uid) = Cell_Cat(thisUnit_f,stem_end_index,'LR');
    %
        UnitsTable.Selectivity_LScene(uid) = Cell_Cat_B(thisUnit_f,stem_end_index,'LScene');
        UnitsTable.Selectivity_RScene(uid) = Cell_Cat_B(thisUnit_f,stem_end_index,'RScene');
        UnitsTable.Selectivity_LR(uid) = Cell_Cat_B(thisUnit_f,stem_end_index,'LR');

%     t=[]; w=[];
%     for f=1:size(thisUnit_f,1)
%         UnitFID = thisUnit_f.ID{f};
%         try
%         [t(f,:),w(f,:)] = CalRDI_Ttest(UnitFID,ROOT,Behav,Spike);
%         catch
%             t=[nan nan nan]; w=[nan nan nan];
%         end
% 
%     end
% 
% 
%     UnitsTable.Selectivity_LScene(uid) = Cell_Cat_C(thisUnit_f,stem_end_index,'LScene',t,w);
%     UnitsTable.Selectivity_RScene(uid) = Cell_Cat_C(thisUnit_f,stem_end_index,'RScene',t,w);
%     UnitsTable.Selectivity_LR(uid) = Cell_Cat_C(thisUnit_f,stem_end_index,'LR',t,w);


    disp([UnitID ' is finished!'])
end

writetable(UnitsTable,[ROOT.Units2 '\UnitsTable_' thisRegion '_forAnalysis.xlsx'],'writemode','replacefile')

sum(UnitsTable.Selectivity_LScene>0 & UnitsTable.Selectivity_LScene<4)
sum(UnitsTable.Selectivity_RScene>0 & UnitsTable.Selectivity_RScene<4)
sum(UnitsTable.Selectivity_LR>0 & UnitsTable.Selectivity_LR<4)

figure
subplot(1,3,1)
histogram(UnitsTable.Selectivity_LScene,'binwidth',0.5)
subplot(1,3,2)
histogram(UnitsTable.Selectivity_RScene,'binwidth',0.5)

subplot(1,3,3)
histogram(UnitsTable.Selectivity_LR,'binwidth',0.5)


%%
function idx = Cell_Cat(units,Dv,p)

if strcmp(p,'LR')
    units(units.PeakBin>Dv,:)=[];
end

if size(units,1)>1
    switch sum(units.(['RateP_' p])<0.05)
        case 0
            idx=0; % all fields are non-selective
        case 1
            idx=2; % one of multiple fields is selective
        otherwise
            units(units.(['RateP_' p])>=0.05,:)=[];
            if prod(units.(['RDI_' p])>0) | prod(units.(['RDI_' p])<0)
                idx=3; % some multiple fields are selective, and homogeneous
            else
                idx=4; % some multiple fields are selective, and heterogeneous
            end
    end
else
    if sum(units.(['RateP_' p])<0.05)
        idx=1; % single field is selective
    else
        idx=0; % single field is non-selective
    end

end
end

%%
function idx = Cell_Cat_B(units,Dv,p)

if strcmp(p,'LR')
    units(units.PeakBin>Dv,:)=[];
end

if size(units,1)>1
    switch sum(abs(units.(['RDI_' p]))>0.1)
        case 0
            idx=0.5; % all fields are non-selective
        case 1
            idx=2; % one of multiple fields is selective
        otherwise
            units(abs(units.(['RDI_' p]))<=0.1,:)=[];
            if prod(units.(['RDI_' p])>0) | prod(units.(['RDI_' p])<0)
                idx=3; % some multiple fields are selective, and homogeneous
            else
                idx=4; % some multiple fields are selective, and heterogeneous
            end
    end
else
    if sum(abs(units.(['RDI_' p]))>0.1)
        idx=1; % single field is selective
    else
        idx=0; % single field is non-selective
    end

end
end
%%
function idx = Cell_Cat_C(units,Dv,p,t,w)



if strcmp(p,'LR')
    units(units.PeakBin>Dv,:)=[];
    c=3;
elseif strcmp(p,'LScene'), c=1;
else, c=2;
end

if size(units,1)>1
    switch sum(t(:,c)<0.05)
        case 0
            idx=0; % all fields are non-selective
        case 1
            idx=2; % one of multiple fields is selective
        otherwise
            units(t(:,c)>=0.05,:)=[];
            if prod(units.(['RDI_' p])>0) | prod(units.(['RDI_' p])<0)
                idx=3; % some multiple fields are selective, and homogeneous
            else
                idx=4; % some multiple fields are selective, and heterogeneous
            end
    end
else
    if sum(t(:,c)<0.05)
        idx=1; % single field is selective
    else
        idx=0; % single field is non-selective
    end

end
end

