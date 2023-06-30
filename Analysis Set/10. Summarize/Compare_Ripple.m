Initial_SWRFilter_common;
warning off
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R3'];
ROOT.Unit = [ROOT.Mother '\Processed Data\units_mat\U0'];
ROOT.Save = [ROOT.Mother '\Processed Data'];
Rip_CA1 = readtable([ROOT.Save '\RipplesTable_' 'CA1_field' '_RDIs_UV_cell_HeteroIn_AllPopul.xlsx']);
Rip_SUB = readtable([ROOT.Save '\RipplesTable_' 'SUB_field' '_RDIs_UV_cell_HeteroIn_AllPopul.xlsx']);
Rip_SUB = readtable([ROOT.Save '\RipplesTable_' 'SUB_refCA1' '_forAnalysis.xlsx']);


pt = table;

Rip_CA1.SpikePerCell = Rip_CA1.spike./Rip_CA1.ensemble;
Rip_SUB.SpikePerCell = Rip_SUB.spike./Rip_SUB.ensemble;


%%
vrList = {'RippleDuration','MeanRaw','MeanFreq','RipplePower','spike','ensemble','SpikePerCell','nFields'};
for v=1:numel(vrList)
vr  = vrList{v};

pt.(vr)(1) = nanmean(Rip_SUB.(vr));
pt.(vr)(2) = nanstd(Rip_SUB.(vr))/sqrt(sum(~isnan(Rip_SUB.(vr))));

pt.(vr)(3) = nanmean(Rip_CA1.(vr));
pt.(vr)(4) = nanstd(Rip_CA1.(vr))/sqrt(sum(~isnan(Rip_CA1.(vr))));

x1 = Rip_SUB.(vr); x1(isnan(x1))=[];
x2 = Rip_CA1.(vr); x2(isnan(x2))=[];
[h,p] = ttest2(x1,x2);

pt.(vr)(5) = p;
end

%%
pt.RipplePower = pt.RipplePower/(10^-12)