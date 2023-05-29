thisRegion2 = 'CA1_field';
RipplesTable_p = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_RDIs_UV' '.xlsx']);
x2 = double(RipplesTable_p.DecodingP_all<0.05);
x1 = [double(RipplesTable_p.pRDI_L_UV<0.05), double(RipplesTable_p.pRDI_R_UV<0.05), double(RipplesTable_p.pRDI_C_UV<0.05)...
    double(min([RipplesTable_p.pRDI_L_UV,RipplesTable_p.pRDI_R_UV,RipplesTable_p.pRDI_C_UV],[],2)<0.05)];

x = x1+2*x2;


edges = [0 0.5 1.5 2.5 3.5];
c = [histcounts(x(:,1),edges);histcounts(x(:,2),edges);histcounts(x(:,3),edges); histcounts(x(:,4),edges)];
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