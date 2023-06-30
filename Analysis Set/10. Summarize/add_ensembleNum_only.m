Initial_SWRFilter_common;
warning off
ROOT.Old = [ROOT.Mother '\Processed Data\ripples_mat\R1'];
ROOT.Save = [ROOT.Mother '\Processed Data\ripples_mat\R2'];
ROOT.Save3 = [ROOT.Mother '\Processed Data\ripples_mat\R3'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];
if ~exist(ROOT.Save), mkdir(ROOT.Save); end
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U2'];

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

thisRegion = 'SUB_refCA1';
RippleList_old = readtable([ROOT.Old '\RipplesTable_' thisRegion '.xlsx'],'ReadRowNames',false);


UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' 'SUB' '_forAnalysis.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' 'SUB_field' '_forAnalysis.xlsx']);

TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

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


%%        %%

for sid=1:size(RippleList_old,1)
    try
        thisRip = RippleList_old(sid,:);
        %                 if thisRip.correctness==1, continue; end
        thisRID = jmnum2str(thisRip.rat,3);
        thisSID = jmnum2str(thisRip.session,2);
        %         if ~strcmp(thisRID,'561'), continue; end

        thisRSID = [jmnum2str(thisRip.rat,3) '-' jmnum2str(thisRip.session,2)];


        if ~strcmp(thisRSID, thisRSID_old)
                  clusters_A = UnitsTable_A(UnitsTable_A.rat==thisRip.rat & UnitsTable_A.session==thisRip.session,:);
     clusters_B = UnitsTable_B(UnitsTable_B.rat==thisRip.rat & UnitsTable_B.session==thisRip.session,:);

                  Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);

            disp([thisRSID ' plotting...'])
            thisRSID_old = thisRSID;

        end

        Clist=jet(256);
        %% Load units

        [aspks_epoch,aspks_epoch_in,aspks_epoch_u,aspks_epoch_u_in,Unitsa,UnitsAa] = LoadUnits(thisRip,clusters_A,Spike,Params,mar);
    [spks_epoch,spks_epoch_in,spks_epoch_u,spks_epoch_u_in,Units,UnitsA] = LoadUnits(thisRip,clusters_B,Spike,Params,mar);

        %%
        RippleList_old.spike(sid) = size(aspks_epoch_in,1);
        RippleList_old.ensemble(sid) = size(aspks_epoch_u_in,1);
        RippleList_old.nFields(sid) = size(spks_epoch_u_in,1);

    end
end
% RippleList = RippleList_old(RippleList_old.ensemble>=3 & strcmp(RippleList_old.experimenter,'LSM'),:)
% RippleList = RippleList(RippleList.nFields>=5,:)
writetable(RippleList_old,[ROOT.Old '\RipplesTable_' thisRegion  '.xlsx'],'WriteMode', 'replacefile')
%%
function  [spks_epoch,spks_epoch_in,spks_epoch_u,spks_epoch_u_in,Units,UnitsA] = LoadUnits(thisRip,clusters,Spike,Params,mar)
spks_epoch=[]; spks_epoch_in=[];
u=0; Units={};UnitsA={};
if ismember('onMazeMaxFR_field', clusters.Properties.VariableNames)
    col4 = 'onMazeMaxFR_field';
else
    col4 = 'onMazeMaxFR';
end

cls_all = size(clusters,1);
for un = 1:cls_all
    [thisRID,thisSID,thisTTID,thisCLID, thisFLID] = parsing_clusterID(clusters.ID{un},1);

    Spk = Spike.(['TT' (thisTTID)]).(['Unit' thisCLID]);

    thisSpks = Spk.t_spk(Spk.t_spk>=thisRip.STtime-mar/Params.Fs & Spk.t_spk<=thisRip.EDtime+mar/Params.Fs);
    thisSpks_in = Spk.t_spk(Spk.t_spk>=thisRip.STtime & Spk.t_spk<=thisRip.EDtime);
    if ~isempty(thisSpks_in)
        
        s1 = ones(size(thisSpks,1),1); s2 = ones(size(thisSpks_in,1),1);
        %                 spks_epoch = [spks_epoch;[thisSpks,s1*un,s1*u,s1*clusters_B.SI(un),...
        %                     s1*clusters_B.RDI_LScene(un), s1*clusters_B.RDI_RScene(un), s1*clusters_B.RDI_LR(un), s1*clusters_A.(['RDI_hetero_' filter_ns])(h)]];
        spks_epoch = [spks_epoch;[thisSpks,s1*un,s1*u,s1*clusters.(col4)(un),...
            s1*clusters.RDI_LScene(un), s1*clusters.RDI_RScene(un), s1*clusters.RDI_LR(un),...
            s1*clusters.Selectivity_LScene(un),s1*clusters.Selectivity_RScene(un),s1*clusters.Selectivity_LR(un)]];
        spks_epoch_in = [spks_epoch_in;[thisSpks_in,s2*un,s2*u,s2*clusters.(col4)(un),...
            s2*clusters.RDI_LScene(un), s2*clusters.RDI_RScene(un), s2*clusters.RDI_LR(un),...
            s2*clusters.Selectivity_LScene(un),s2*clusters.Selectivity_RScene(un),s2*clusters.Selectivity_LR(un)]];
        % spks_epoch_in = [spks_epoch_in;[thisSpks_in,s2*un,s2*u,s2*clusters_B.SI(un),...
        %                     s2*clusters_B.RDI_LScene(un), s2*clusters_B.RDI_RScene(un), s2*clusters_B.RDI_LR(un),s2*clusters_A.(['RDI_hetero_' filter_ns])(h)]];
        u=u+1;
        if ~isnan(str2double(thisFLID))
        Units = [Units; [thisTTID '-' thisCLID '-' thisFLID]];
        else
             Units = [Units; [thisTTID '-' thisCLID]];
        end
        UnitsA = [UnitsA; [clusters.ID(un)]];
    end

end
if ~isempty(spks_epoch_in)
[~,ia,~] = unique(spks_epoch(:,3),'rows');
spks_epoch_u = spks_epoch(ia,:);
[~,ia,~] = unique(spks_epoch_in(:,3),'rows');
spks_epoch_u_in = spks_epoch_in(ia,:);
else
    spks_epoch_u=[];
    spks_epoch_u_in=[];
end
end
