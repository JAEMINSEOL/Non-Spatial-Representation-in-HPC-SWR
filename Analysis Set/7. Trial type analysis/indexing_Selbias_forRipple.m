Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Save '\ripples_mat\R4_10ms'];
ROOT.Fig3 = [ROOT.Save '\ripples_mat\ProfilingSheet\R11_AI_hot'];
ROOT.Units = [ROOT.Save '\units_mat\U1'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];

dir = '';
ROOT.Fig = [ROOT.Fig3];
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


thisRegion = 'CA1';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable = UnitsTable_A;
TT_table = readtable([ROOT.Info '\TT_table.xlsx']);

Experimenter = {'LSM','JS','SEB'};
unit = 5.0000e-04; ti = 1;

RipplesTable_all = table;
thisRSID_old = '';
mar = 0.05*Params.Fs;
dur = 0.4*Params.Fs;
thisFRMapSCALE=2;
Params.tbinDuration = 0.005;


for sid=1:size(RipplesTable,1)

    if RipplesTable.pRDI_L(sid)<0.05
        if RipplesTable.mRDI_L(sid)>0
            RipplesTable.SceneBias(sid) = 1;
            RipplesTable.Decoding_CxtMatch(sid) = Decoding_match(RipplesTable.DecodingP_Z(sid));
        else
            RipplesTable.SceneBias(sid) = 3;
            RipplesTable.Decoding_CxtMatch(sid) = Decoding_match(RipplesTable.DecodingP_B(sid));
        end

    elseif RipplesTable.pRDI_R(sid)<0.05
         if RipplesTable.mRDI_R(sid)>0
            RipplesTable.SceneBias(sid) = 2;
            RipplesTable.Decoding_CxtMatch(sid) = Decoding_match(RipplesTable.DecodingP_P(sid));
        else
            RipplesTable.SceneBias(sid) = 4;
            RipplesTable.Decoding_CxtMatch(sid) = Decoding_match(RipplesTable.DecodingP_M(sid));
         end
    else
        RipplesTable.SceneBias(sid) = 0;
    end

    if RipplesTable.pRDI_C(sid)<0.05
        if RipplesTable.mRDI_C(sid)>0
            RipplesTable.ChoiceBias(sid)=1;
        else
            RipplesTable.ChoiceBias(sid)=2;
        end
    else
        RipplesTable.ChoiceBias(sid)=0;
    end
end
RipplesTable.Decoding_All = RipplesTable.DecodingP_all<0.05;


 writetable(RipplesTable,[ROOT.Save '\RipplesTable_' thisRegion '_forAnalysis_RDI.xlsx']);
%%
ParsedRipTable = RipplesTable(RipplesTable.ChoiceBias>0,:);

figure;
count = hist2d([ParsedRipTable.ChoiceBias,mod(ParsedRipTable.context,2)+1],[0:2],[0:2]);
for i=1:size(count,1)
    count(i,:) = count(i,:)/sum(count(i,:));
end
imagesc(count)
colormap jet

ParsedRipTable = RipplesTable(RipplesTable.SceneBias>0,:);
figure;
count = hist2d([ParsedRipTable.SceneBias,ParsedRipTable.context],[0:4],[0:4]);
for i=1:size(count,2)
    count(:,i) = count(:,i)/sum(count(:,i));
end
imagesc(count)
colormap jet
title('Replay')

figure;
count = hist2d([ParsedRipTable.SceneBias,ParsedRipTable.Decoding_CxtMatch],[0:4],[0:1]);
imagesc(count)

function p = Decoding_match(decod)
if decod<0.05
    p=1;
else
    p=0;
end
end