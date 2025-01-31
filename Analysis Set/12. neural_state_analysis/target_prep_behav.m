addpath(genpath('D:\HPC-SWR project\Analysis Program'))


Initial_SWRFilter_common;
warning off
ROOT.Processed = [ROOT.Mother '\Processed Data'];
ROOT.Save = [ROOT.Processed];
ROOT.Behav = [ROOT.Processed '\behavior_mat'];
%%
clear RipplesTable

Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);

RegionList = {'CA1','SUB'};
Target=struct;
binsz=0.2;

%
window=7;
alpha=1.75;
smooth=1;
FR_Thre=0.1;
variance=80;
mode = 2;
gaussFilt = gausswin(window,alpha);
clear gaussFilter
colorMap = parula(25600);

for i=1:(window+1)/2
    gaussFilter(:,i)= gaussFilt / sum(gaussFilt((window+1)/2+1-i:window));
end
%% add unit ratio
% for sid=1:size(SessionList,1)
for sid=82:86

    thisRID = SessionList.rat(sid);
    thisSID = SessionList.session(sid);
    thisRSID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];
    %         load([ROOT.Behav '\' thisRSID '.mat']);
    BehavTable = readtable([ROOT.Behav '\' thisRSID '.xlsx']);
    if ~exist([ROOT.Behav '\' thisRSID '.mat']), continue; end
    load([ROOT.Behav '\' thisRSID '.mat']);
    Behav.t2 = Behav.t-Behav.t(1);
    Behav.v2 = Behav.velocity.speed(knnsearch(Behav.t,Behav.t2));
    for reg=1:2
        try

            thisR = RegionList{reg};
            thisRegion = thisR;
            TargetTT = GetTargetTT(ROOT,thisRSID,thisRegion,Params,1);
            Spike=LoadSpikeData(ROOT, thisRSID, [1:24], Params.cellfindn);
            targ =[];

            bins = linspace(Behav.t2(1),Behav.t2(end),round(max(Behav.t2)/binsz)+1);

            bin_pos=[];
            for b = 1:length(bins)-1
                bin_pos(b,1) = nanmedian(Behav.y_linearized(Behav.t2>=bins(b) & Behav.t2<bins(b+1)));
            end
            %         bin_time=[];
            %         for b = 1:length(bins)-1
            %             bin_time(b) = sum(Behav.y>=bins(b) & Behav.y<bins(b+1))/30;
            %         end

            for t=1:size(TargetTT,1)
                thisTT = TargetTT(t);
                if isfield(Spike,['TT' num2str(thisTT)])
                    field_list = fieldnames(Spike.(['TT' num2str(thisTT)]));
                    for fid = 1:size(field_list,1)
                        thisUnit = str2double(field_list{fid}(5:end));
                        try
                            thisspks = Spike.(['TT' num2str(thisTT)]).(field_list{fid}).t_spk-Behav.t(1);
                            thisMap = load([ROOT.Raw.Map '\rat' thisRSID '-' num2str(TargetTT(t)) '-' jmnum2str(thisUnit,2)],'skaggsMap1D');
                            thisMap = thisMap.skaggsMap1D{1};

                            bin_spk=[];
                            for b = 1:length(bins)-1
                                bin_spk(b) = sum(thisspks>=bins(b) & thisspks<bins(b+1));
                            end
                            %         targ = [targ,thisMap];
                            targ = [targ,bin_spk'./0.5];
                        end
                    end
                end
            end



            id = isnan(targ(:,1)) | bin_pos>108;
            targ(id,:)=[];
            targ(:,isnan(targ(1,:)))=[];
            targ = zscore(targ);
            c=bin_pos;
            c(id,:)=[];
            iso = isomap(targ,3);
            gp = gplvm(targ,3);

            Target.(['rat' jmnum2str(thisRID,3) '_' jmnum2str(thisSID,2) '_' thisR]) = targ;
            Isomap.(['rat' jmnum2str(thisRID,3) '_' jmnum2str(thisSID,2) '_' thisR]) = iso;
            GPLVM.(['rat' jmnum2str(thisRID,3) '_' jmnum2str(thisSID,2) '_' thisR]) = gp;
            Pos.(['rat' jmnum2str(thisRID,3) '_' jmnum2str(thisSID,2) '_' thisR]) = c;
        end
    end

end
save([ROOT.Save '\processed_behav2.mat'],'Target','Isomap','GPLVM','Pos')

%%
load([ROOT.Save '\processed_behav2.mat'],'Target','Isomap','GPLVM','Pos')
%%
ROOT.Fig = [ROOT.Save '\neural_state\isomap_behav'];
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end
for sid=64:67
    try
        reg=2;
        thisRID = SessionList.rat(sid);
        thisSID = SessionList.session(sid);
        thisRSID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];
        % thisRSID = '561-02';
        thisRID = thisRSID(1:3); thisSID=thisRSID(end-1:end); thisR = RegionList{reg};

        targ=Target.(['rat' thisRID '_' thisSID '_' thisR]);
        gp = GetPC_GaussianFiltering(Isomap.(['rat' thisRID '_' thisSID '_' thisR]), window, gaussFilter);
        pos = Pos.(['rat' thisRID '_' thisSID '_' thisR]);
        % pos=c;
        figure('position',[0 0 1000 1000])
        scatter3(gp(:,1),gp(:,2),gp(:,3),20,pos,'filled')
        colormap jet
        colorbar
        xlabel('PC1'); ylabel('PC2'); zlabel('PC3')
        % xlim([-20 20]); xlim([-20 20])
        title([thisRSID ',' thisR ', Isomap, behav'])
        saveas(gca,[ROOT.Fig '\' thisRSID '.fig'])
        saveas(gca,[ROOT.Fig '\' thisRSID '.png'])
    end
end
% close all
%%
ROOT.Fig = [ROOT.Save '\neural_state\isomap_behav'];
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end
for sid=64:88

    for reg=1:2
        try
            % thisRSID = '415-13'; reg=2;
            % thisRID = thisRSID(1:3); thisSID=thisRSID(end-1:end); thisR = RegionList{reg};
            thisRID = SessionList.rat(sid);
            thisSID = SessionList.session(sid);
            thisRSID = [jmnum2str(thisRID,3) '-' jmnum2str(thisSID,2)];
            % thisRSID = '561-02';
            thisRID = thisRSID(1:3); thisSID=thisRSID(end-1:end); thisR = RegionList{reg};

            load([ROOT.Behav '\' thisRSID '.mat']);
            Behav.t2 = Behav.t-Behav.t(1);
            Behav.v2 = Behav.velocity.speed(knnsearch(Behav.t,Behav.t2));

            bins = linspace(Behav.t2(1),Behav.t2(end),round(max(Behav.t2)/binsz)+1);
            bin_pos=[];
            for b = 1:length(bins)-1
                bin_pos(b,1) = nanmedian(Behav.y_linearized(Behav.t2>=bins(b) & Behav.t2<bins(b+1)));
            end

            v=[20 50 40];
            targ=Target.(['rat' thisRID '_' thisSID '_' thisR]);
            iso=GetPC_GaussianFiltering(Isomap.(['rat' thisRID '_' thisSID '_' thisR]), window, gaussFilter);
            pos = Pos.(['rat' thisRID '_' thisSID '_' thisR]);
            idx = bin_pos>108;

            id = knnsearch(Behav.t2,bins');
            c = Behav.cont(id,:); c(Behav.correctness(:,2),:)=1;
            c=c(~idx,:);
            % iso=gp;

            figure('position',[0 0 1000 1000]); hold on

            ax1=axes;
            id = c(:,2);
            scatter3(ax1,iso(id,1),iso(id,2),iso(id,3),20,pos(id,:),'filled')
            % ax1.Visible = 'off';
            colormap(ax1,'winter')
            view(ax1,v)
            xlabel('PC1'); ylabel('PC2'); zlabel('PC3')

            ax2=axes;
            id = c(:,4);
            scatter3(ax2,iso(id,1),iso(id,2),iso(id,3),20,pos(id,:),'filled')
            ax2.Visible = 'off'; colormap(ax2,'autumn')
            view(ax2,v)
            linkaxes([ax1,ax2])
            set([ax1,ax2],'Position',[.17 .11 .685 .815]);
            cb1 = colorbar(ax1,'Position',[.05 .11 .03 .815]);
            cb2 = colorbar(ax2,'Position',[.9 .11 .03 .815]);


            % figure; hold on
            % id = c(:,1);
            % scatter3(iso(id,1),iso(id,2),iso(id,3),20,'b','filled')
            % id = c(:,2);
            % scatter3(iso(id,1),iso(id,2),iso(id,3),20,'r','filled')


            sgtitle([thisRSID ',' thisR ', Pebbles vs. Mountain, Isomap'])

            saveas(gca,[ROOT.Fig '\' thisRSID '_right.fig'])
            saveas(gca,[ROOT.Fig '\' thisRSID '_right.png'])
        end
    end
end