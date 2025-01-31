clear RipplesTable

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

RegionList = {'CA1','SUB'};
SessionTable=struct;
%%
for reg=1:2

thisR = RegionList{reg};
thisRegion0 = thisR;
thisRegion = thisR;
thisRegion2 = [thisR '_field'];

SessionTable.(thisRegion0)=table;

R1 = readtable([ROOT.Save '\RipplesTable_' thisRegion0 '_final.xlsx']);
U1 = readtable([ROOT.Save '\UnitsTable_' thisRegion '_final.xlsx']);

R1_sleep = readtable([ROOT.Save '_sleep\RipplesTable_' thisRegion0 '_final.xlsx']);
U1_sleep = readtable([ROOT.Save '_sleep\UnitsTable_' thisRegion '_final.xlsx']);

R0 = readtable([ROOT.Mother '\Processed Data\ripples_mat\R2\RipplesTable_Behav_' thisRegion0 '_speed_filtered.xlsx']);
R0_sleep = readtable([ROOT.Mother '\Processed Data_sleep\ripples_mat\R2\RipplesTable_Behav_' thisRegion0 '_filtered.xlsx']);

B0 = readtable([ROOT.Mother '\Processed Data_sleep\sleep_epoch.xlsx']);


%%
for sid = 1:size(SessionList)
    thisRID = SessionList.rat(sid);
    thisSID = SessionList.session(sid);
    thisRSID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];

    thisU1 = U1(U1.rat==thisRID & U1.session==thisSID,:);

    thisR1 = R1(R1.rat==thisRID & R1.session==thisSID,:);
    thisR0 = R0(R0.rat==thisRID & R0.session==thisSID & ~(R0.Overlap==0),:);

        thisR1_sleep = R1_sleep(R1_sleep.rat==thisRID & R1_sleep.session==thisSID,:);
    thisR0_sleep = R0_sleep(R0_sleep.rat==thisRID & R0_sleep.session==thisSID,:);
    dist = thisB1.duration_ITI;

    dist(dist>mean(dist)+std(dist))=[];

    if ~isempty(thisR0)
        temp=table;
            thisB1 = readtable([ROOT.Mother '\Processed Data\behavior_mat\' thisRSID '.xlsx']);
    thisB0 = B0(strcmp(B0.RSID,thisRSID),:);
    temp.ID{1}= thisRSID;
    temp.numRips_behav(1) = size(thisR0,1);
    temp.dur_behav = sum(dist);
    temp.RipRate_behav = temp.numRips_behav/temp.dur_behav;

    temp.numReacts_behav = size(thisR1,1);
    temp.ReactRate_behav = temp.numReacts_behav/temp.dur_behav;

%     temp.numRips_pre = sum(thisR0_sleep.sleep_pre_post==1);
%     temp.dur_pre = thisB0.pre_end-thisB0.pre_start;
%     temp.RipRate_pre = temp.numRips_pre/temp.dur_pre;
% 
%     temp.numRips_post = sum(thisR0_sleep.sleep_pre_post==2);
%     temp.dur_post = thisB0.post_end-thisB0.post_start;
%     temp.RipRate_post = temp.numRips_post/temp.dur_post;

    temp.numUnits = size(thisU1,1);

    SessionTable.(thisRegion0) = [SessionTable.(thisRegion0);temp];
    end

end
end

%%
writetable(SessionTable.SUB,[ROOT.Save '\SessionTable_SUB.xlsx'],'writemode','replacefile')
writetable(SessionTable.CA1,[ROOT.Save '\SessionTable_CA1.xlsx'],'writemode','replacefile')
