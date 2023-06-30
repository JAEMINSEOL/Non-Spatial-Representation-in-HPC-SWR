thisRegion2 = 'CA1_field';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_RDIs_UV_cell_HeteroIn_AllPopul' '.xlsx']);


% RipplesTable=RipplesTable(RipplesTable.correctness==1,:);
% RipplesTable=RipplesTable(RipplesTable.RippleDuration>0.03,:);
% 
%  writetable(RipplesTable,[ROOT.Save '\RipplesTable_' thisRegion2 '_RDIs_UV_cell_HeteroIn' '.xlsx'],'writemode','replacefile')

x2 = double(RipplesTable.DecodingP_all<0.05);
x1 = [double(RipplesTable.pRDI_L_UV<0.05), double(RipplesTable.pRDI_R_UV<0.05), ...
    double(min([RipplesTable.pRDI_L_UV,RipplesTable.pRDI_R_UV],[],2)<0.05),...
    double(RipplesTable.pRDI_C_UV<0.05)...
    double(min([RipplesTable.pRDI_L_UV,RipplesTable.pRDI_R_UV,RipplesTable.pRDI_C_UV],[],2)<0.05)];

x = x1+2*x2;


edges = [0 0.5 1.5 2.5 3.5];
c = [histcounts(x(:,1),edges);histcounts(x(:,2),edges);histcounts(x(:,3),edges); histcounts(x(:,4),edges); histcounts(x(:,5),edges)];
%%
[n,e] = histcounts(x(:,1),edges)
figure;
subplot(1,4,1)
pie(c(1,[2 1 3 4]))
subplot(1,4,2)
pie(c(2,[2 1 3 4]))
subplot(1,4,3)
pie(c(3,[2 1 3 4]))

subplot(1,4,4)
pie(c(4,[2 1 3 4]))

%%
sum(abs(units.RDI_LScene)>0.1 | abs(units.RDI_RScene)>0.1)
sum(abs(units.RDI_LR)>0.1)

nsid = nanmin([RipplesTable.pRDI_L_UV,RipplesTable.pRDI_R_UV,RipplesTable.pRDI_C_UV],[],2)<0.05;
sum(nsid)

histogram(abs(UnitsTable_B.RDI_LR),'binwidth',0.025)
histogram(UnitsTable_A.Selectivity_LR,'binwidth',0.5)