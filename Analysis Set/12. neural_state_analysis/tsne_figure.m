
cd('D:\HPC-SWR project\Processed Data')
addpath(genpath('D:\HPC-SWR project\Analysis Program'))


Initial_SWRFilter_common;
warning off
ROOT.Processed = [ROOT.Mother '\Processed Data'];
ROOT.Save = [ROOT.Processed '\temp_final'];
alg = 'TSNE';
gradname = 'TDA_Ratio';

ROOT.Fig = [ROOT.Processed '\neural_state\' alg '_' gradname];
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end

pca_l = load([ROOT.Save '\processed_L.mat']);
pca_r = load([ROOT.Save '\processed_R.mat']);
pca_c = load([ROOT.Save '\processed_C.mat']);
%%
window=7;
alpha=1.75;
smooth=1;
FR_Thre=0.1;
variance=80;
mode = 2;
gaussFilt = gausswin(window,alpha);
clear gaussFilter
colorMap = parula(25600);

sz=20;
[grad,im]=colorGradient(hex2rgb('f4a026'),hex2rgb('4a5ca8'),256);

for i=1:(window+1)/2
    gaussFilter(:,i)= gaussFilt / sum(gaussFilt((window+1)/2+1-i:window));
end
%% 3D scatter_sample
% CList = {'999999','ff8192','80a0ff','30bf72','2f4858'};
% field_name = 'r415_s13_SUB';
% 
% 
% targ = [zscore(Target.(field_name))];
% 
% % initial_dims = 3; % t-SNE의 차원 (예: 2D)
% % initial_data = pca(targ', 'NumComponents', initial_dims);
% thisPCs = tsne(targ,'NumDimensions',3,'InitialY',initial_data,'InitialY',PC.(field_name)(:,1:3));
% time = Rips.(field_name).STtime; time=time-time(1);
% for rid=1:size(Rips.(field_name),1)
%     Rips.(field_name).trial_num(rid) = str2double(Rips.(field_name).trial{rid}(end-2:end));
% end
% cls = Rips.(field_name).StartTime_fromTrialEnd;
% if Dimension.(field_name) > 2
%     fig=figure('position',[266,100,800,800]);
%     hold on;
% 
%     idx_z = find(Rips.(field_name).context==1);
%     idx_p = find(Rips.(field_name).context==2);
%     idx_b = find(Rips.(field_name).context==3);
%     idx_m = find(Rips.(field_name).context==4);
% 
%     scatter3(thisPCs(idx_z,3),thisPCs(idx_z,2),thisPCs(idx_z,1),sz*1.5,hex2rgb(CList{2}),'filled');
%     scatter3(thisPCs(idx_p,3),thisPCs(idx_p,2),thisPCs(idx_p,1),sz*1.5,hex2rgb(CList{3}),'filled');
%     scatter3(thisPCs(idx_b,3),thisPCs(idx_b,2),thisPCs(idx_b,1),sz*1.5,hex2rgb(CList{4}),'filled');
%     scatter3(thisPCs(idx_m,3),thisPCs(idx_m,2),thisPCs(idx_m,1),sz*1.5,hex2rgb(CList{5}),'filled');
% 
%     cluster_id = Rips.(field_name).context;
%     cvi1 = silhouette_index(thisPCs, cluster_id);
%     grid on
%     g=gca; g.FontSize=12;
%     xlabel('PC3'); ylabel('PC2'); zlabel('PC1');
%     title([field_name ', cvi= ' num2str(cvi1)],'Interpreter','none');
%     view([40,10,10])
% 
% end

%% 3D scatter_batch_trial condition
% cvi_table = table;
% CList = {'aaaaaa','F3C623','4FBFD6','E72416','4F63AD','ff8192','80a0ff'};
% 
% Rips = pca_l.Rips;
% Target = pca_l.Target;
% Dimension = pca_l.Dimension;
% field_list = fieldnames(Target);
% 
% for fid = 1:size(field_list,1)
%     try
% 
%         field_name = field_list{fid};
% 
%         if contains(field_name,'sleep_pre') || contains(field_name,'sleep_post'), continue; else, cvi_table.sleep(fid)=0; end
%         if contains(field_name,'SUB'), cvi_table.region(fid)=1; else, cvi_table.region(fid)=2; end
% 
%         
% %         timeTE = Rips.(field_name).StartTime_fromTrialEnd; 
% %         idx = timeTE<10;
% 
%         thisPCs.L = pca_l.(alg).(field_name);
%         thisPCs.R = pca_r.(alg).(field_name);
%         thisPCs.C = pca_c.(alg).(field_name);
% time = Rips.(field_name).STtime; 
%         time=time-time(1);
% %         if Dimension.(field_name) > 2
% cvi_table.nRips(fid) = size(Rips,1);
% 
%             %% trial condition label
% %             fig=figure('position',[-1000,20,3000,800]);
% %             v=[20 40 50];
% %             subplot(1,3,1)
% %             hold on;
% %             scatter3(thisPCs(:,1),thisPCs(:,2),thisPCs(:,3),sz*3,time,'filled')
% %             grid on
% %             g=gca; g.FontSize=12;
% %             xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% %             title(field_name,'Interpreter','none');
% % %             view([40,10,10])
% %             c = colorbar('location','eastoutside');
% %             c.Label.String = 'Time (s)';
% %             colormap jet
% %             view(v)
% % 
% %             subplot(1,3,2)
% %             hold on;
% % 
% %             nfields='_max';
% %             idx_z = find(Rips.(field_name).context==1);
% %             idx_p = find(Rips.(field_name).context==2);
% %             idx_b = find(Rips.(field_name).context==3);
% %             idx_m = find(Rips.(field_name).context==4);
% % 
% % %             scatter3(thisPCs(:,3),thisPCs(:,2),thisPCs(:,1),sz,hex2rgb(CList{1}))
% % 
% % 
% %             scatter3(thisPCs(idx_z,1),thisPCs(idx_z,2),thisPCs(idx_z,3),sz*1.5,hex2rgb(CList{2}),'filled');
% %             scatter3(thisPCs(idx_p,1),thisPCs(idx_p,2),thisPCs(idx_p,3),sz*1.5,hex2rgb(CList{3}),'filled');
% %             scatter3(thisPCs(idx_b,1),thisPCs(idx_b,2),thisPCs(idx_b,3),sz*1.5,hex2rgb(CList{4}),'filled');
% %             scatter3(thisPCs(idx_m,1),thisPCs(idx_m,2),thisPCs(idx_m,3),sz*1.5,hex2rgb(CList{5}),'filled');
% %             view(v)
% % 
% %             mean_z=mean(thisPCs(idx_z,1:3,1));
% %             mean_p=mean(thisPCs(idx_p,1:3,1));
% %             mean_b=mean(thisPCs(idx_b,1:3,1));
% %             mean_m=mean(thisPCs(idx_m,1:3,1));
% % 
% % 
% %             cluster_id = Rips.(field_name).context;
% %             cvi1 = silhouette_index(thisPCs, cluster_id);
% %             cvi3 = dunns(4,thisPCs, cluster_id);
% % 
% %             grid on
% %             g=gca; g.FontSize=12;
% %             xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% %             title([field_name ', cvi= ' num2str(cvi3)],'Interpreter','none');
% % %             view([40,10,10])
% % 
% %             legend({'Zebra','Pebbles','Bamboo','Mountain'})
% % 
% %           subplot(1,3,3)
% %             hold on;
% %             idx_l = find(Rips.(field_name).context==1 | Rips.(field_name).context==3);
% %             idx_r = find(Rips.(field_name).context==2 | Rips.(field_name).context==4);
% % 
% %             % scatter3(thisPCs(:,3),thisPCs(:,2),thisPCs(:,1),sz,hex2rgb(CList{1}))
% %             scatter3(thisPCs(idx_l,1),thisPCs(idx_l,2),thisPCs(idx_l,3),sz*1.5,hex2rgb(CList{2}),'filled');
% %             scatter3(thisPCs(idx_r,1),thisPCs(idx_r,2),thisPCs(idx_r,3),sz*1.5,hex2rgb(CList{3}),'filled');
% %             view(v)
% % 
% %             mean_l=mean(thisPCs(idx_l,1:3,1));
% %             mean_r=mean(thisPCs(idx_r,1:3,1));
% % 
% %             %             scatter3(mean_l(3),mean_l(2),mean_l(1),sz*10,'r','filled','pentagram','markeredgecolor','y')
% %             %             scatter3(mean_r(3),mean_r(2),mean_r(1),sz*10,'b','filled','pentagram','markeredgecolor','y')
% %             %
% %             %
% %             %             plot3([mean_l(3) mean_r(3)],[mean_l(2) mean_r(2)],[mean_l(1) mean_r(1)],'k:','LineWidth',2)
% % 
% %             %             d2 = pdist([mean_l;mean_r]);
% % 
% %             cluster_id = mod(Rips.(field_name).context,2);
% %             cvi2 = silhouette_index(thisPCs, cluster_id);
% %              cvi4 = dunns(2,thisPCs, cluster_id);
% % 
% %             grid on
% %             g=gca; g.FontSize=12;
% %             xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% %             title([field_name ', cvi= ' num2str(cvi4)],'Interpreter','none');
% % %             view([40,10,10])
% % 
% %             legend({'Left','Right'})
% % 
% %             cvi_table.silhS(fid) = cvi1;
% %             cvi_table.silhC(fid) = cvi2;
% % 
% %                         cvi_table.dunnS(fid) = cvi3;
% %             cvi_table.dunnC(fid) = cvi4;
% % cvi_table.hopkins(fid) = hopkins(thisPCs,10);
% %% TDA label
%             fig=figure('position',[266,244,3000,800]);
%             subplot(1,3,1)
%             hold on;
%             [nb,nz,mb,mz] = color_scatter_2('L',Rips,field_name,{CList{1} CList{2} CList{4}},sz,thisPCs.L);
%             legend({'None','Bamboo','Zebra'})
% 
%             subplot(1,3,2)
%             hold on;
%             [nm,np,mm,mp] = color_scatter_2('R',Rips,field_name,{CList{1} CList{3} CList{5}},sz,thisPCs.R);
%             legend({'None','Mountain','Pebbles'})
% 
%             subplot(1,3,3)
%             hold on;
%             [nl,nr,ml,mr] = color_scatter_2('C',Rips,field_name,{CList{1} CList{6} CList{7}},sz,thisPCs.C);
%             legend({'None','Left','Right'})
% 
%             %             scatter3(mean_z(3),mean_z(2),mean_z(1),sz*10,'r','filled','pentagram','markeredgecolor','y')
%             %             scatter3(mean_p(3),mean_p(2),mean_p(1),sz*10,'b','filled','pentagram','markeredgecolor','y')
%             %             scatter3(mean_b(3),mean_b(2),mean_b(1),sz*10,hex2rgb('5f9750'),'filled','pentagram','markeredgecolor','y')
%             %             scatter3(mean_m(3),mean_m(2),mean_m(1),sz*10,'k','filled','pentagram','markeredgecolor','y')
%             %
%             %
%             %             plot3([mean_z(3) mean_p(3)],[mean_z(2) mean_p(2)],[mean_z(1) mean_p(1)],'k:','LineWidth',2)
%             %             plot3([mean_z(3) mean_b(3)],[mean_z(2) mean_b(2)],[mean_z(1) mean_b(1)],'k:','LineWidth',2)
%             %             plot3([mean_z(3) mean_m(3)],[mean_z(2) mean_m(2)],[mean_z(1) mean_m(1)],'k:','LineWidth',2)
%             %             plot3([mean_b(3) mean_p(3)],[mean_b(2) mean_p(2)],[mean_b(1) mean_p(1)],'k:','LineWidth',2)
%             %             plot3([mean_m(3) mean_p(3)],[mean_m(2) mean_p(2)],[mean_m(1) mean_p(1)],'k:','LineWidth',2)
%             %             plot3([mean_b(3) mean_m(3)],[mean_b(2) mean_m(2)],[mean_b(1) mean_m(1)],'k:','LineWidth',2)
% 
%             %             d1 = mean([pdist([mean_z;mean_p]),pdist([mean_z;mean_b]),pdist([mean_z;mean_m]),pdist([mean_b;mean_p]),pdist([mean_m;mean_p]),pdist([mean_b;mean_m])]);
% 
%             %     text(m(1),m(2),m(3),['dist= ' num2str(pdist([mean_pos;mean_neg]))])
% 
%             sgtitle([field_name ',' alg],'Interpreter','none')
%   
% %%
%             saveas(gca,[ROOT.Fig '\' field_name '.fig'])
%             saveas(gca,[ROOT.Fig '\' field_name '.png'])
% 
% %         end
% 
% 
%         close all
%     catch
%         disp([field_name ' failed'])
%         close all
% 
%     end
% 
% end

%% 3D scatter_batch_gradient
cvi_table = table;
CList = {'aaaaaa','F3C623','4FBFD6','E72416','4F63AD','ff8192','80a0ff'};

ctrial=[];
for c=1:4
ctrial(c,:) = hex2rgb(CList{c+1});
end

Rips = pca_l.Rips;
Target = pca_l.Target;
Dimension = pca_l.Dimension;
field_list = fieldnames(Target);

for fid = 1:size(field_list,1)
        try

        field_name = field_list{fid};

        if contains(field_name,'sleep_pre') || contains(field_name,'sleep_post'), continue; else, cvi_table.sleep(fid)=0; end
        if contains(field_name,'SUB'), cvi_table.region(fid)=1; else, cvi_table.region(fid)=2; end

%                     fig = openfig([ROOT.Fig '\' field_name '.fig']);
%             fig.Position = [266,244,3000,800];
%             ax = findall(gcf, 'Type', 'axes');
%             view(ax(1),[40,10,10])
%             view(ax(2),[40,10,10])
%             view(ax(3),[40,10,10])
            
        timeTE = Rips.(field_name).StartTime_fromTrialEnd; 
        idx = timeTE<10;

        thisPCs.L = pca_l.(alg).(field_name);
        thisPCs.R = pca_r.(alg).(field_name);
        thisPCs.C = pca_c.(alg).(field_name);
time = Rips.(field_name).STtime; 
        time=time-time(1);
%         if Dimension.(field_name) > 2
cvi_table.nRips(fid) = size(Rips.(field_name),1);

        
% TDA label
            fig=figure('position',[266,244,3000,800]);
            subplot('position',[.05 .1 .27 .8])
            hold on;
            [km,kn] = color_scatter_3('L',Rips,field_name,sz,thisPCs.L);
            legend({'None-TDA','Zebra-TDA','Bamboo-TDA'})
            title(['for left scene pair, knn = ' jjnum2str(kn,3) ', kmeans = ' jjnum2str(km,3)])
            cvi_table.kn(fid,1) = kn; cvi_table.kmean(fid,1) = km;
            %             colormap(grad)
%             colormap(ctrial)

            subplot('position',[.37 .1 .27 .8])
            hold on;
            [km,kn] = color_scatter_3('R',Rips,field_name,sz,thisPCs.R);
            legend({'None-TDA','Pebbles-TDA','Mountain-TDA'})
            title(['for right scene pair, knn = ' jjnum2str(kn,3) ', kmeans = ' jjnum2str(km,3)])
            cvi_table.kn(fid,2) = kn; cvi_table.kmean(fid,2) = km;
            %             colormap(grad)
%             colormap(ctrial)

            subplot('position',[.69 .1 .25 .8])
            hold on;
            [km,kn] = color_scatter_3('C',Rips,field_name,sz,thisPCs.C);
            legend({'None-TDA','Left-TDA','Right-TDA'})
            title(['for choice pair, knn = ' jjnum2str(kn,3) ', kmeans = ' jjnum2str(km,3)])
            cvi_table.kn(fid,3) = kn; cvi_table.kmean(fid,3) = km;
            %             colormap(grad)
%             colormap(ctrial)


            %             scatter3(mean_z(3),mean_z(2),mean_z(1),sz*10,'r','filled','pentagram','markeredgecolor','y')
            %             scatter3(mean_p(3),mean_p(2),mean_p(1),sz*10,'b','filled','pentagram','markeredgecolor','y')
            %             scatter3(mean_b(3),mean_b(2),mean_b(1),sz*10,hex2rgb('5f9750'),'filled','pentagram','markeredgecolor','y')
            %             scatter3(mean_m(3),mean_m(2),mean_m(1),sz*10,'k','filled','pentagram','markeredgecolor','y')
            %
            %
            %             plot3([mean_z(3) mean_p(3)],[mean_z(2) mean_p(2)],[mean_z(1) mean_p(1)],'k:','LineWidth',2)
            %             plot3([mean_z(3) mean_b(3)],[mean_z(2) mean_b(2)],[mean_z(1) mean_b(1)],'k:','LineWidth',2)
            %             plot3([mean_z(3) mean_m(3)],[mean_z(2) mean_m(2)],[mean_z(1) mean_m(1)],'k:','LineWidth',2)
            %             plot3([mean_b(3) mean_p(3)],[mean_b(2) mean_p(2)],[mean_b(1) mean_p(1)],'k:','LineWidth',2)
            %             plot3([mean_m(3) mean_p(3)],[mean_m(2) mean_p(2)],[mean_m(1) mean_p(1)],'k:','LineWidth',2)
            %             plot3([mean_b(3) mean_m(3)],[mean_b(2) mean_m(2)],[mean_b(1) mean_m(1)],'k:','LineWidth',2)

            %             d1 = mean([pdist([mean_z;mean_p]),pdist([mean_z;mean_b]),pdist([mean_z;mean_m]),pdist([mean_b;mean_p]),pdist([mean_m;mean_p]),pdist([mean_b;mean_m])]);

            %     text(m(1),m(2),m(3),['dist= ' num2str(pdist([mean_pos;mean_neg]))])

            sgtitle([field_name ',' alg ',' gradname],'Interpreter','none')
  
%%
            saveas(gca,[ROOT.Fig '\' field_name '.fig'])
            saveas(gca,[ROOT.Fig '\' field_name '.svg'])
            saveas(gca,[ROOT.Fig '\' field_name '.png'])
            exportgraphics(gcf,[ROOT.Fig '\' field_name '.pdf'],'BackgroundColor','none','ContentType','vector')

%         end


        close all
    catch
        disp([field_name ' failed'])
        close all

    end

end

%%
writetable(cvi_table,[ROOT.Processed '\neural_state\CVI_table_' alg '_' gradname '.xlsx'],'writemode','replacefile');
%%
cvi_table = readtable([ROOT.Processed '\neural_state\CVI_table_' alg '_' gradname '.xlsx']);
cvi_temp = cvi_table(cvi_table.sleep==0 & cvi_table.nRips>100,:);


figure
subplot(1,2,2); hold on

dat1 = reshape([cvi_temp.kn_1(cvi_temp.region==1,:);cvi_temp.kn_2(cvi_temp.region==1,:);cvi_temp.kn_3(cvi_temp.region==1,:)],1,[]);
dat2 = reshape([cvi_temp.kn_1(cvi_temp.region==2,:);cvi_temp.kn_2(cvi_temp.region==2,:);cvi_temp.kn_3(cvi_temp.region==2,:)],1,[]);

% m1 = mean(dat1);
% s1 = std(dat1)/sqrt(numel(dat1));
% m2 = mean(dat2);
% s2 = std(dat2)/sqrt(numel(dat2));
% bar([1 2],[m1 m2])
% errorbar([m1 m2],[s1 s2],'color','k')

boxplot([dat1 dat2],[ones(length(dat1),1); ones(length(dat2),1)*2])
scatter(ones(length(dat1),1),dat1,'filled')
scatter(ones(length(dat2),1)*2,dat2,'filled')
[p,h,stats] = ranksum(dat1,dat2);
title(['p=' num2str(p) ',Z=' num2str(stats.zval)])

xticks([1 2]); xticklabels({'SUB','CA1'})
ylabel('K-NN cluster quality')
ylim([0 1])
g=gca; g.FontSize=12;

subplot(1,2,1); hold on

dat3 = reshape([cvi_temp.kmean_1(cvi_temp.region==1,:);cvi_temp.kmean_2(cvi_temp.region==1,:);cvi_temp.kmean_3(cvi_temp.region==1,:)],1,[]);
dat4 = reshape([cvi_temp.kmean_1(cvi_temp.region==2,:);cvi_temp.kmean_2(cvi_temp.region==2,:);cvi_temp.kmean_3(cvi_temp.region==2,:)],1,[]);
% m3 = mean(dat3);
% s1 = std(dat3)/sqrt(numel(dat3));
% m4 = mean(dat4);
% s2 = std(dat4)/sqrt(numel(dat4));
% bar([1 2],[m3 m4])
% errorbar([m3 m4],[s1 s2],'color','k')

boxplot([dat3 dat4],[ones(length(dat3),1); ones(length(dat4),1)*2])
scatter(ones(length(dat3),1),dat3,'filled')
scatter(ones(length(dat4),1)*2,dat4,'filled')
[p,h,stats] = ranksum(dat3,dat4);
title(['p=' num2str(p) ',Z=' num2str(stats.zval)])

xticks([1 2]); xticklabels({'SUB','CA1'})
ylabel('K-means cluster quality')
ylim([0 1])
g=gca; g.FontSize=12;

sgtitle([alg ',' gradname],'Interpreter','none')
% saveas(gca,[ROOT.Processed '\neural_state\' 'cluster_qual_' alg '_' gradname '.png'])
% saveas(gca,[ROOT.Processed '\neural_state\' 'cluster_qual_' alg '_' gradname '.svg'])
%%
% cvi_temp = cvi_table(cvi_table.hopkins~=0 & cvi_table.nRips>40,:);
% 
% 
% figure
% subplot(1,1,1); hold on
% 
% m1 = mean(cvi_temp.hopkins(cvi_temp.region==1));
% s1 = std(cvi_temp.hopkins(cvi_temp.region==1))/sqrt(sum(cvi_temp.region==1));
% m2 = mean(cvi_temp.hopkins(cvi_temp.region==2));
% s2 = std(cvi_temp.hopkins(cvi_temp.region==2))/sqrt(sum(cvi_temp.region==2));
% bar([1 2],[m1 m2])
% errorbar([m1 m2],[s1 s2],'color','k')
% 
% xticks([1 2]); xticklabels({'SUB','CA1'})
% ylabel('neural state difference (Scene)')
% % ylim([-.05 .05])
% g=gca; g.FontSize=12;
% 
% 
% xticks([1 2]); xticklabels({'SUB','CA1'})
% ylabel('hopkins statistics')
% % ylim([-.05 .05])
% g=gca; g.FontSize=12;
% %%
% t1 = cvi_temp.hopkins(cvi_temp.region==1);
% t2 = cvi_temp.hopkins(cvi_temp.region==2);
% 
% [h,p] = ranksum(t1,t2)
% 
% %%
% % 데이터 준비
% X = rand(100, 3); % 예제 데이터
% X=  thisPCs;
% 
% % Hopkins 통계량 계산
% m = 10; % 샘플 크기
% H_value = hopkins(X, m);
% fprintf('Hopkins Statistic: %.4f\n', H_value);

figure
scatter(dat1,dat3,'filled');
hold on
scatter(dat2,dat4,'filled')

xlabel('K-NN cluster quality'); ylabel('K-means cluster quality')
legend('SUB','CA1')
%%
%% function
function [n_pos,n_neg,mean_pos,mean_neg] = color_scatter_2(ch,Rips,field_name,CList,sz,thisPCs)
nfields='_max';
thisRips=Rips.(field_name);
idx_p = find(thisRips.(['pBinomDev_' ch '_UV'])<0.05 & thisRips.(['nRDI_' ch nfields])>=5 & thisRips.(['Ratio_' ch'])>.5);
idx_n = find(thisRips.(['pBinomDev_' ch '_UV'])<0.05 & thisRips.(['nRDI_' ch nfields])>=5 & thisRips.(['Ratio_' ch'])<.5);

scatter3(thisPCs(:,3),thisPCs(:,2),thisPCs(:,1),sz,hex2rgb(CList{1}))
scatter3(thisPCs(idx_p,3),thisPCs(idx_p,2),thisPCs(idx_p,1),sz*1.5,hex2rgb(CList{2}),'filled');
scatter3(thisPCs(idx_n,3),thisPCs(idx_n,2),thisPCs(idx_n,1),sz*1.5,hex2rgb(CList{3}),'filled');

n_pos = size(idx_p,1);
n_neg = size(idx_n,1);
try
    try
        if n_pos>1
            mean_pos=mean(thisPCs(idx_p,1:3,1));
% AllPoints = thisPCs(idx_p,1:3,1);
%             K = 1;
% for idx=1:length(AllPoints)
% [~,r] = findNearestNeighbors(AllPoints,AllPoints(idx,:),K);
% density(idx) = 1/(4*pi*r.^3/3);
% end

        else
            mean_pos=(thisPCs(idx_p,1:3,1));
        end

        scatter3(mean_pos(3),mean_pos(2),mean_pos(1),sz*10,hex2rgb(CList{2}),'filled','pentagram','markeredgecolor','k')
    end
    try

        if n_neg>1
            mean_neg=mean(thisPCs(idx_n,1:3,1));
        else
            mean_neg=(thisPCs(idx_n,1:3,1));
        end
        scatter3(mean_neg(3),mean_neg(2),mean_neg(1),sz*10,hex2rgb(CList{3}),'filled','pentagram','markeredgecolor','k')
    end

    plot3([mean_pos(3) mean_neg(3)],[mean_pos(2) mean_neg(2)],[mean_pos(1) mean_neg(1)],'k:','LineWidth',2)

    m=mean([mean_pos(3) mean_neg(3);mean_pos(2) mean_neg(2);mean_pos(1) mean_neg(1)],2);
    d = pdist([mean_pos;mean_neg]);

    %     text(m(1),m(2),m(3),['dist= ' num2str(pdist([mean_pos;mean_neg]))])
catch
    d=0;
end

grid on
g=gca; g.FontSize=12;
xlabel('PC3'); ylabel('PC2'); zlabel('PC1');
title([ 'dist= ' num2str(d)],'Interpreter','none');
view([40,10,10])
end

function [km,kn]=color_scatter_3(ch,Rips,field_name,sz,thisPCs)
nfields='_max';
thisRips=Rips.(field_name);
idx_p = find(thisRips.(['pBinomDev_' ch '_UV'])<0.05 & thisRips.(['nRDI_' ch nfields])>=5 & thisRips.(['Ratio_' ch'])>.5);
idx_n = find(thisRips.(['pBinomDev_' ch '_UV'])<0.05 & thisRips.(['nRDI_' ch nfields])>=5 & thisRips.(['Ratio_' ch'])<.5);

clab = thisRips.(['Ratio_' ch]);
% clab = thisRips.(['StartTime_fromTrialEnd'])/15;
% clab = thisRips.(['context']);
% clab = thisRips.(['pBinomDev_' ch '_UV']);

scatter3(thisPCs(:,3),thisPCs(:,2),thisPCs(:,1),sz*1.5,clab,'filled')
% scatter3(thisPCs(idx_n,3),thisPCs(idx_n,2),thisPCs(idx_n,1),sz*7,'b','diamond','filled','markeredgecolor','k');
% scatter3(thisPCs(idx_p,3),thisPCs(idx_p,2),thisPCs(idx_p,1),sz*7,'r','diamond','filled','markeredgecolor','k');
% scatter3(thisPCs(idx_n,3),thisPCs(idx_n,2),thisPCs(idx_n,1),sz*7,thisRips.(['pBinomDev_' ch '_UV'])(idx_n),'diamond','filled','markeredgecolor','k');
% scatter3(thisPCs(idx_p,3),thisPCs(idx_p,2),thisPCs(idx_p,1),sz*7,thisRips.(['pBinomDev_' ch '_UV'])(idx_p),'pentagram','filled','markeredgecolor','k');
% scatter3(thisPCs(idx_n,3),thisPCs(idx_n,2),thisPCs(idx_n,1),sz*7,thisRips.context(idx_n),'diamond','filled','markeredgecolor','k');
% scatter3(thisPCs(idx_p,3),thisPCs(idx_p,2),thisPCs(idx_p,1),sz*7,thisRips.context(idx_p),'pentagram','filled','markeredgecolor','k');

colormap(jet)
clim([0 1])
cb = colorbar('location','eastoutside');

% cb.Ticks = [1:1:4];
% cb.TickLabels = {'Zebra','Pebbles','Bamboo','Mountains'};

cb.Ticks = [0 .5 1];
cb.TickLabels = {'-1', '0', '1'};

grid on
g=gca; g.FontSize=12;
xlabel('PC3'); ylabel('PC2'); zlabel('PC1');
view([40,10,10])

kni=[]; kmi=[]; acci=[];
k = 2;
data = thisPCs;
labels = clab;
tole = .05;

for i=1:1000
% KNN 모델 생성
rid = randperm(size(data,1),round(size(data,1)/5));
knnModel = fitcknn(data(rid,:), labels(rid));

% 클러스터링 결과 예측
predictedLabels = predict(knnModel, data);


% 실제 레이블과 예측된 레이블 비교
kni(i) = sum(predictedLabels == labels) / length(labels);
acci(i) = mean(abs((labels-predictedLabels)) <=tole);

% K-means 클러스터링 수행
[idx,~] = kmeans(data, k);

% 실루엣 계수 계산
silhouetteValues = silhouette(data, idx);
kmi(i) = mean(silhouetteValues);
end

km=mean(kmi);
kn=mean(kni);
kn = mean(acci);


end
%%
