%% UnitSpec_compare
Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Processed ''];
ROOT.Rip0 = [ROOT.Save '\ripples_mat\R0'];
ROOT.Rip = [ROOT.Save '\ripples_mat\R2'];
ROOT.Rip4 = [ROOT.Save '\ripples_mat\R4_SUB_refCA1'];
ROOT.Rip5 = [ROOT.Save '\ripples_mat\R5_cell_AllPopul'];
ROOT.Fig = [ROOT.Save '\ripples_mat\ProfilingSheet\R25_ca1'];
ROOT.Unit1 = [ROOT.Save '\units_mat\U1'];
ROOT.Units = [ROOT.Save '\units_mat\U2'];
ROOT.Behav = [ROOT.Save '\behavior_mat'];


CList = [ [207 8 23]/255;[23 84 181]/255];


RegionList = {'SUB','CA1'};

load([ROOT.Processed '\Manuscript\ripple_SI_4.mat']);
% RipplesTable_SI.CA1 = Rip_CA1;

Rip_SUB = RipplesTable_SI.SUB;
Rip_CA1 = RipplesTable_SI.CA1;

% Rip_temp2 = readtable([ROOT.Save '\RipplesTable_' 'CA1' '_forAnalysis_final.xlsx']);
% save([ROOT.Processed '\Manuscript\ripple_SI.mat'],'RipplesTable_SI');
%% for fig 4
crit_list = [0.05:0.05:0.4];
unit_list = {'SUB','CA1'};

figure('Position',[86,165,1215,813])
for ulid = 1:2
    subplot(1,2,ulid); hold on
for cri = 1:length(crit_list)
    crit_si =crit_list(cri);
    crit_binom=0.05;
    Rip_temp = RipplesTable_SI.(unit_list{ulid}).(['SI_' jmnum2str(crit_si*100,3)]);

    Rip_temp.pBinomDev_L_UV(Rip_temp.nRDI_L_max<5) = nan;
    Rip_temp.pBinomDev_R_UV(Rip_temp.nRDI_R_max<5) = nan;
    Rip_temp.pBinomDev_C_UV(Rip_temp.nRDI_C_max<5) = nan;
    Rip_temp(Rip_temp.rat==232 & Rip_temp.session==4,:) = [];

    sum_all = sum((Rip_temp.pBinomDev_L_UV<crit_binom | Rip_temp.pBinomDev_R_UV<crit_binom)  | (Rip_temp.pBinomDev_C_UV<crit_binom ));
%     sum_all = sum(max([Rip_temp.nRDI_L_max Rip_temp.nRDI_R_max Rip_temp.nRDI_C_max],[],2)>=5);
    sum_ssi(cri) = sum((Rip_temp.pBinomDev_L_UV<crit_binom | Rip_temp.pBinomDev_R_UV<crit_binom)  & ~(Rip_temp.pBinomDev_C_UV<crit_binom ))/sum_all;
    sum_csi(cri) = sum(~(Rip_temp.pBinomDev_L_UV<crit_binom | Rip_temp.pBinomDev_R_UV<crit_binom)  & (Rip_temp.pBinomDev_C_UV<crit_binom ))/sum_all;
    sum_cssi(cri) = sum((Rip_temp.pBinomDev_L_UV<crit_binom | Rip_temp.pBinomDev_R_UV<crit_binom)  & (Rip_temp.pBinomDev_C_UV<crit_binom ))/sum_all;
end

plot(sum_ssi,'color',hex2rgb('e8aa42'),'Marker','o')
plot(sum_csi,'color',hex2rgb('d83f31'),'Marker','o')
plot(sum_cssi,'color',hex2rgb('025464'),'Marker','o')

xticks([1:length(crit_list)])
xticklabels(crit_list)
xlabel('Selective index criterion of cells')
ylim([0 1])
ylabel(' proportion in Task-related SWRs')

legend({'Scene only','Choice only','Scene and Choice'},'location','northoutside','Orientation','horizontal')

title(unit_list{ulid})
end

%% for fig 4
crit_list = [0.05:0.05:0.4];
unit_list = {'SUB','CA1'};

figure('Position',[86,165,1215,813])
for ulid = 1:2
    subplot(1,2,ulid); hold on
for cri = 1:length(crit_list)
    crit_si =crit_list(cri);
    crit_binom=0.05;
    Rip_temp = RipplesTable_SI.(unit_list{ulid}).(['SI_' jmnum2str(crit_si*100,3)]);

    Rip_temp.pBinomDev_L_UV(Rip_temp.nRDI_L_max<5) = nan;
    Rip_temp.pBinomDev_R_UV(Rip_temp.nRDI_R_max<5) = nan;
    Rip_temp.pBinomDev_C_UV(Rip_temp.nRDI_C_max<5) = nan;
    Rip_temp(Rip_temp.rat==232 & Rip_temp.session==4,:) = [];


    sum_all = sum(max([Rip_temp.nRDI_L_max Rip_temp.nRDI_R_max Rip_temp.nRDI_C_max],[],2)>=5);
%     sum_all = sum(Rip_temp.nPCs>=5);
%     sum_all = sum(Rip_temp.nPCs>0);
    sum_ssi(cri) = sum((Rip_temp.pBinomDev_L_UV<crit_binom | Rip_temp.pBinomDev_R_UV<crit_binom))/sum_all;
    sum_csi(cri) = sum((Rip_temp.pBinomDev_C_UV<crit_binom ))/sum_all;
    sum_cssi(cri) = sum((Rip_temp.pBinomDev_C_UV<crit_binom | Rip_temp.pBinomDev_L_UV<crit_binom | Rip_temp.pBinomDev_R_UV<crit_binom ))/sum_all;
end

plot(sum_ssi,'color',hex2rgb('e8aa42'),'Marker','o')
plot(sum_csi,'color',hex2rgb('d83f31'),'Marker','o')
plot(sum_cssi,'color',hex2rgb('000000'),'Marker','o')

xticks([1:length(crit_list)])
xticklabels(crit_list)
xlabel('Selective index criterion of cells')
ylim([0 1])
ylabel('Task-related SWR proportion')

legend({'Scene','Choice','Scene or Choice'},'location','northoutside','Orientation','horizontal')

title(unit_list{ulid})
end

    
