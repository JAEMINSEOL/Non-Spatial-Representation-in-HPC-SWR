addpath(genpath('D:\HPC-SWR project\Analysis Program'))


Initial_SWRFilter_common;
warning off
ROOT.Processed = [ROOT.Mother '\Processed Data'];
ROOT.Save = [ROOT.Processed];

ROOT.Fig = [ROOT.Save '\neural_state\trial_condition'];
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end

load([ROOT.Save '\processed_pca_fr.mat'])
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

for i=1:(window+1)/2
    gaussFilter(:,i)= gaussFilt / sum(gaussFilt((window+1)/2+1-i:window));
end
%% 3D scatter_sample
% field_name = 'r415_s13_SUB';
% 
% thisPCs = PC.([field_name '_GaussFiltered']);
% time = Rips.(field_name).STtime; time=time-time(1);
% for rid=1:size(Rips.(field_name),1)
%     Rips.(field_name).trial_num(rid) = str2double(Rips.(field_name).trial{rid}(end-2:end));
% end
% cls = Rips.(field_name).StartTime_fromTrialEnd;
% if Dimension.(field_name) > 2
%     fig=figure('position',[266,100,800,800]);
%     hold on;
%     scatter3(thisPCs(:,3),thisPCs(:,2),thisPCs(:,1),sz,cls,'filled')
%     grid on
%     g=gca; g.FontSize=12;
%     xlabel('PC3'); ylabel('PC2'); zlabel('PC1');
%     title(field_name,'Interpreter','none');
%     view([40,10,10])
% 
% end

%% 3D scatter_batch

com_table=table; dist_table = table;
% %%%%%%%%%%%%%%%%%%%%%
CList = {'999999','ff8192','80a0ff'};

field_name = 'r415_s13_SUB';
field_list = fieldnames(Target);

for fid = 1:size(field_list,1)
    try

        field_name = field_list{fid};
        if contains(field_name,'sleep_pre'), com_table.sleep(fid)=1; elseif contains(field_name,'sleep_post'), com_table.sleep(fid)=2; else, com_table.sleep(fid)=0; end
        if contains(field_name,'SUB'), com_table.region(fid)=1; else, com_table.region(fid)=2; end

        if contains(field_name,'sleep_pre'), dist_table.sleep(fid)=1; elseif contains(field_name,'sleep_post'), dist_table.sleep(fid)=2; else, dist_table.sleep(fid)=0; end
        if contains(field_name,'SUB'), dist_table.region(fid)=1; else, dist_table.region(fid)=2; end

        thisPCs = PC.([field_name '_GaussFiltered']);
        time = Rips.(field_name).STtime; time=time-time(1);
        if Dimension.(field_name) > 2
            fig=figure('position',[266,244,800,800]);
            hold on;
            scatter3(thisPCs(:,3),thisPCs(:,2),thisPCs(:,1),sz*3,time,'filled')
            grid on
            g=gca; g.FontSize=12;
            xlabel('PC3'); ylabel('PC2'); zlabel('PC1');
            title(field_name,'Interpreter','none');
            view([40,10,10])
            c = colorbar('location','eastoutside');
            c.Label.String = 'Time (s)';
            saveas(gca,[ROOT.Save '\' field_name '_time.fig'])
            saveas(gca,[ROOT.Save '\neural_state\' field_name '_time.png'])


            fig=figure('position',[266,244,3000,800]);
            subplot(1,3,1)
            hold on;
            [nb,nz,mb,mz] = color_scatter_2('L',Rips,field_name,CList,sz,thisPCs);
            legend({'None','Bamboo','Zebra'})

            subplot(1,3,2)
            hold on;
            [nm,np,mm,mp] = color_scatter_2('R',Rips,field_name,CList,sz,thisPCs);
            legend({'None','Mountain','Pebbles'})

            subplot(1,3,3)
            hold on;
            [nl,nr,ml,mr] = color_scatter_2('C',Rips,field_name,CList,sz,thisPCs);
            legend({'None','Left','Right'})


            com_table.n_Bamboo(fid) = nb;
            com_table.n_Zebra(fid) = nz;
            com_table.n_Mountain(fid) = nm;
            com_table.n_Pebbles(fid) = np;
            com_table.n_Left(fid) = nl;
            com_table.n_Right(fid) = nr;

                        com_table.COM_Bamboo{fid} = mb;
            com_table.COM_Zebra{fid} = mz;
            com_table.COM_Mountain{fid} = mm;
            com_table.COM_Pebbles{fid} = mp;
            com_table.COM_Left{fid} = ml;
            com_table.COM_Right{fid} = mr;

            if isempty(pdist([mb;mz])) || min([nb,nz])<2, dist_table.distL(fid) = nan; else, dist_table.distL(fid) = pdist([mb;mz]); end
            if isempty(pdist([mm;mp])) || min([nm,np])<2, dist_table.distR(fid) = nan; else, dist_table.distR(fid) = pdist([mm;mp]); end
            if isempty(pdist([ml;mr])) || min([nl,nr])<2, dist_table.distC(fid) = nan; else, dist_table.distC(fid) = pdist([ml;mr]); end

            saveas(gca,[ROOT.Fig '\' field_name '_TDA.fig'])
            saveas(gca,[ROOT.Fig '\' field_name '_TDA.png'])


        end

        

        close all
    catch
        disp([field_name ' failed'])
        close all

    end

end

%% 3D scatter_batch_trial condition
dist_table = table;
% %%%%%%%%%%%%%%%%%%%%%
CList = {'999999','ff8192','80a0ff','30bf72','2f4858'};

field_name = 'r415_s13_SUB';
field_list = fieldnames(Target);

for fid = 1:size(field_list,1)
    try

        field_name = field_list{fid};
        if contains(field_name,'sleep_pre'), com_table.sleep(fid)=1; elseif contains(field_name,'sleep_post'), com_table.sleep(fid)=2; else, com_table.sleep(fid)=0; end
        if contains(field_name,'SUB'), com_table.region(fid)=1; else, com_table.region(fid)=2; end

        if contains(field_name,'sleep_pre'), dist_table.sleep(fid)=1; elseif contains(field_name,'sleep_post'), dist_table.sleep(fid)=2; else, dist_table.sleep(fid)=0; end
        if contains(field_name,'SUB'), dist_table.region(fid)=1; else, dist_table.region(fid)=2; end

        thisPCs = PC.([field_name '_GaussFiltered']);
        time = Rips.(field_name).STtime; time=time-time(1);
        if Dimension.(field_name) > 2

            fig=figure('position',[266,244,3000,800]);

subplot(1,3,1)
            hold on;
            scatter3(thisPCs(:,3),thisPCs(:,2),thisPCs(:,1),sz*3,time,'filled')
            grid on
            g=gca; g.FontSize=12;
            xlabel('PC3'); ylabel('PC2'); zlabel('PC1');
            title(field_name,'Interpreter','none');
            view([40,10,10])
            c = colorbar('location','eastoutside');
            c.Label.String = 'Time (s)';


            subplot(1,3,2)
            hold on;

            nfields='_max';
            idx_z = find(Rips.(field_name).context==1);
            idx_p = find(Rips.(field_name).context==2);
            idx_b = find(Rips.(field_name).context==3);
            idx_m = find(Rips.(field_name).context==4);

            % scatter3(thisPCs(:,3),thisPCs(:,2),thisPCs(:,1),sz,hex2rgb(CList{1}))
            scatter3(thisPCs(idx_z,3),thisPCs(idx_z,2),thisPCs(idx_z,1),sz*1.5,hex2rgb(CList{2}),'filled');
            scatter3(thisPCs(idx_p,3),thisPCs(idx_p,2),thisPCs(idx_p,1),sz*1.5,hex2rgb(CList{3}),'filled');
            scatter3(thisPCs(idx_b,3),thisPCs(idx_b,2),thisPCs(idx_b,1),sz*1.5,hex2rgb(CList{4}),'filled');
            scatter3(thisPCs(idx_m,3),thisPCs(idx_m,2),thisPCs(idx_m,1),sz*1.5,hex2rgb(CList{5}),'filled');


            mean_z=mean(thisPCs(idx_z,1:3,1));
            mean_p=mean(thisPCs(idx_p,1:3,1));
            mean_b=mean(thisPCs(idx_b,1:3,1));
            mean_m=mean(thisPCs(idx_m,1:3,1));

            scatter3(mean_z(3),mean_z(2),mean_z(1),sz*10,'r','filled','pentagram','markeredgecolor','y')
            scatter3(mean_p(3),mean_p(2),mean_p(1),sz*10,'b','filled','pentagram','markeredgecolor','y')
            scatter3(mean_b(3),mean_b(2),mean_b(1),sz*10,hex2rgb('5f9750'),'filled','pentagram','markeredgecolor','y')
            scatter3(mean_m(3),mean_m(2),mean_m(1),sz*10,'k','filled','pentagram','markeredgecolor','y')


            plot3([mean_z(3) mean_p(3)],[mean_z(2) mean_p(2)],[mean_z(1) mean_p(1)],'k:','LineWidth',2)
            plot3([mean_z(3) mean_b(3)],[mean_z(2) mean_b(2)],[mean_z(1) mean_b(1)],'k:','LineWidth',2)
            plot3([mean_z(3) mean_m(3)],[mean_z(2) mean_m(2)],[mean_z(1) mean_m(1)],'k:','LineWidth',2)
            plot3([mean_b(3) mean_p(3)],[mean_b(2) mean_p(2)],[mean_b(1) mean_p(1)],'k:','LineWidth',2)
            plot3([mean_m(3) mean_p(3)],[mean_m(2) mean_p(2)],[mean_m(1) mean_p(1)],'k:','LineWidth',2)
            plot3([mean_b(3) mean_m(3)],[mean_b(2) mean_m(2)],[mean_b(1) mean_m(1)],'k:','LineWidth',2)

            d1 = mean([pdist([mean_z;mean_p]),pdist([mean_z;mean_b]),pdist([mean_z;mean_m]),pdist([mean_b;mean_p]),pdist([mean_m;mean_p]),pdist([mean_b;mean_m])]);

            %     text(m(1),m(2),m(3),['dist= ' num2str(pdist([mean_pos;mean_neg]))])


            grid on
            g=gca; g.FontSize=12;
            xlabel('PC3'); ylabel('PC2'); zlabel('PC1');
            title([field_name ', dist= ' num2str(d1)],'Interpreter','none');
            view([40,10,10])

            legend({'Zebra','Pebbles','Bamboo','Mountain'})


            subplot(1,3,3)
            hold on;
            idx_l = find(Rips.(field_name).context==1 | Rips.(field_name).context==3);
            idx_r = find(Rips.(field_name).context==2 | Rips.(field_name).context==4);

            % scatter3(thisPCs(:,3),thisPCs(:,2),thisPCs(:,1),sz,hex2rgb(CList{1}))
            scatter3(thisPCs(idx_l,3),thisPCs(idx_l,2),thisPCs(idx_l,1),sz*1.5,hex2rgb(CList{2}),'filled');
            scatter3(thisPCs(idx_r,3),thisPCs(idx_r,2),thisPCs(idx_r,1),sz*1.5,hex2rgb(CList{3}),'filled');

            mean_l=mean(thisPCs(idx_l,1:3,1));
            mean_r=mean(thisPCs(idx_r,1:3,1));

            scatter3(mean_l(3),mean_l(2),mean_l(1),sz*10,'r','filled','pentagram','markeredgecolor','y')
            scatter3(mean_r(3),mean_r(2),mean_r(1),sz*10,'b','filled','pentagram','markeredgecolor','y')


            plot3([mean_l(3) mean_r(3)],[mean_l(2) mean_r(2)],[mean_l(1) mean_r(1)],'k:','LineWidth',2)

            d2 = pdist([mean_l;mean_r]);

            grid on
            g=gca; g.FontSize=12;
            xlabel('PC3'); ylabel('PC2'); zlabel('PC1');
            title([field_name ', dist= ' num2str(d2)],'Interpreter','none');
            view([40,10,10])

            legend({'Left','Right'})

            dist_table.distS(fid) = d1;
           dist_table.distC(fid) = d2;

            saveas(gca,[ROOT.Fig '\' field_name '_TDA.fig'])
            saveas(gca,[ROOT.Fig '\' field_name '_TDA.png'])

        end

        
        close all
    catch
        disp([field_name ' failed'])
        close all

    end

end
%%
dist_table.mean_dist = nanmean([dist_table.distL,dist_table.distR,dist_table.distC],2);

dist_temp = dist_table(dist_table.sleep==0,:);

figure
boxplot(dist_temp.mean_dist,dist_temp.region)
hold on
scatter(dist_temp.region,dist_temp.mean_dist)

xticklabels({'SUB','CA1'})
ylabel('mean neural state difference')
%%
sub = com_table.distL(com_table.sleep==0 & com_table.region==1);
sub(sub==0)=[];
%%
% 
% field_name = 'r415_s13_SUB';
% 
% if contains(field_name,'sleep'), nfields='_max'; else, nfields='_UV'; end
% 
% Zebra.([field_name '_GaussFiltered'])=mean(PC.([field_name '_GaussFiltered'])(idx_p,1:Dimension.(field_name)),1);
% 
% 
% 
% for i=1:size(PC.([field_name '_GaussFiltered']),1)
%     x=PC.([field_name '_GaussFiltered'])(i,1:Dimension.(field_name));
%     y=Zebra.([field_name '_GaussFiltered']);
%     z=Variance.(field_name)(1:Dimension.(field_name));
%     EUdist.([field_name '_Weighted_GaussFiltered'])(i,:)=GetWeightedEuclideanDistance(x,y,z);
%     EUdist.([field_name '_GaussFiltered'])(i,:)=pdist2(x,y);
% end

%% function
function [n_pos,n_neg,mean_pos,mean_neg] = color_scatter_2(ch,Rips,field_name,CList,sz,thisPCs)
nfields='_max';
idx_p = find(Rips.(field_name).(['pBinomDev_' ch '_UV'])<0.05 & Rips.(field_name).(['nRDI_' ch nfields])>=5 & Rips.(field_name).(['Ratio_' ch'])>.5);
idx_n = find(Rips.(field_name).(['pBinomDev_' ch '_UV'])<0.05 & Rips.(field_name).(['nRDI_' ch nfields])>=5 & Rips.(field_name).(['Ratio_' ch'])<.5);

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

        scatter3(mean_pos(3),mean_pos(2),mean_pos(1),sz*10,'r','filled','pentagram','markeredgecolor','r')
    end
    try

        if n_neg>1
            mean_neg=mean(thisPCs(idx_n,1:3,1));
        else
            mean_neg=(thisPCs(idx_n,1:3,1));
        end
        scatter3(mean_neg(3),mean_neg(2),mean_neg(1),sz*10,'b','filled','pentagram','markeredgecolor','b')
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
title([field_name ', dist= ' num2str(d)],'Interpreter','none');
view([40,10,10])
end

