addpath(genpath('D:\HPC-SWR project\Analysis Program'))


Initial_SWRFilter_common;
warning off
ROOT.Processed = [ROOT.Mother '\Processed Data'];
ROOT.Save = [ROOT.Processed '\temp_final'];

load([ROOT.Save '\processed_pca.mat'])
rng(50)



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
TList={'L','R','C'};
field_list = fieldnames(TargetF);
%%

for ti = 1:3
for fid = 1:size(field_list,1)

field_name = field_list{fid};
thisR = Rips.(field_name); thisU = Units.(field_name); thisF = Fields.(field_name); targ = TargetF.(field_name);
thisF.RDI_L = thisF.RDI_LScene; thisF.RDI_R = thisF.RDI_RScene; thisF.RDI_C = thisF.RDI_LR;
if size(targ,1)<20, continue; end
% for t=1:size(Target.(field_name),1)
%     Target.(field_name)(t,:) = Target.(field_name)(t,:) / Rips.(field_name).RippleDuration(t);
% end

for r=1:size(thisR,1)
    thisR.trial_num(r) = str2double(thisR.trial{r}(end-2:end));
end
% targ = [zscore(Target.(field_name))];
% targ = [targ, thisR.RippleDuration];


targ = double(targ>=1);

for u=1:size(thisU,1)
       for r=1:size(thisR,1)
          idf = find(contains(thisF.ID,thisU.ID{u}));
        if thisR.(['Ratio_' TList{ti}])(r)>=0.5
            [~,ia] = max(thisF.(['RDI_' TList{ti}])(idf));
            idf(ia)=[];
            targ(r,idf)=0;
        elseif thisR.(['Ratio_' TList{ti}])(r)<=0.5
                        [~,ia] = min(thisF.(['RDI_' TList{ti}])(idf));
            idf(ia)=[];
            targ(r,idf)=0;
        end
       end
        if thisU.(['RDI_' TList{ti} '_min'])(u)>-0.1 && thisU.(['RDI_' TList{ti} '_max'])(u)<0.1
        targ(:,idf)=nan;
        end
end

inan = isnan(targ(1,:));
targ(:,inan)=[]; thisF(inan,:)=[];

% raj=rand(size(targ))*10^(-10);
%  targ=targ+raj;

if size(targ,2)<3, continue; end

[coef.(field_name), PC.(field_name),~,~,Variance.(field_name)] = pca(targ);
PC.([field_name '_GaussFiltered']) = GetPC_GaussianFiltering(PC.(field_name), window, gaussFilter);
TSNE.(field_name) = tsne(targ,'NumDimensions',3,'InitialY',PC.(field_name)(:,1:3));
ISOMAP.(field_name) = isomap(GetPC_GaussianFiltering(targ, window, gaussFilter),3);
GPLVM.(field_name) =  gplvm(GetPC_GaussianFiltering(targ, window, gaussFilter),3);


for i=1:length(Variance.(field_name))
    if sum(Variance.(field_name)(1:i)) > variance
        Dimension.(field_name)=i;
        break;
    end
end
% figure;
% thisPCs = ISOMAP.(field_name);
% scatter3(thisPCs(:,3),thisPCs(:,2),thisPCs(:,1),20,'k')
% title(field_name)
end

% cmdscale(zscore(Target.(field_name)(:,:)))

save([ROOT.Save '\processed_' TList{ti} '.mat'],'Target','TargetF','Rips','Units','Fields','coef','PC','TSNE','Variance','Dimension','ISOMAP','GPLVM')
%%
end
% tsne_figure;
load([ROOT.Save '\processed_pca.mat'])
%%
TList={'L','R','C'};
cmap = [1 1 1; 0 0 0; 1 0 0; 0 0 1];
ROOT.Fig = [ROOT.Processed '\neural_state\Input_data'];
for fid = 1:size(field_list,1)

field_name = field_list{fid};
TTList = {'Left','Right','Choice'};
for ti=1:2
thisR = Rips.(field_name); thisU = Units.(field_name); thisF = Fields.(field_name); targ = TargetF.(field_name);
targ = double(targ>=1);
for u=1:size(thisU,1)

     for r=1:size(thisR,1)
         idf = find(contains(thisF.ID,thisU.ID{u}));
         if thisR.(['Ratio_' TList{ti}])(r)>=0.5
             if thisR.(['pBinomDev_' TList{ti} '_UV'])(r)<0.05 & thisR.(['nRDI_' TList{ti} '_max'])(r)>=5
                 targ(r,targ(r,:)==1)=2;
             end
             [~,ia] = max(thisF.(['RDI_' TList{ti} 'Scene'])(idf));
             idf(ia)=[];
             targ(r,idf)=0;
             [~,ia] = min(thisF.(['RDI_' TList{ti} 'Scene'])(idf));
            idf(ia)=[];
            targ(r,idf)=nan;
         elseif thisR.(['Ratio_' TList{ti}])(r)<=0.5
             if thisR.(['pBinomDev_' TList{ti} '_UV'])(r)<0.05 & thisR.(['nRDI_' TList{ti} '_max'])(r)>=5
                 targ(r,targ(r,:)==1)=3;
             end
            [~,ia] = min(thisF.(['RDI_' TList{ti} 'Scene'])(idf));
            idf(ia)=[];
            targ(r,idf)=0;
                        [~,ia] = max(thisF.(['RDI_' TList{ti} 'Scene'])(idf));
            idf(ia)=[];
            targ(r,idf)=nan;
         end
                     
    end
    idf = find(contains(thisF.ID,thisU.ID{u}));
    ia = abs(thisF.(['RDI_' TList{ti} 'Scene'])(idf))<0.1;
    targ(:,idf(ia))=nan;


        if thisU.(['RDI_' TList{ti} '_min'])(u)>-0.1 && thisU.(['RDI_' TList{ti} '_max'])(u)<0.1
        targ(:,idf)=nan;
        end

end

inan = isnan(targ(1,:));
targ(:,inan)=[]; thisF(inan,:)=[];
[r,ia] = sort(thisR.(['Ratio_' TList{ti}]));
[u,ib] = sort(thisF.(['RDI_' TList{ti} 'Scene']));
targ = targ(ia,ib);
figure('position',[115,100,1198,859]);

h=heatmap(targ); grid off
colormap(cmap)
clim([0 3])
title(field_name)
h.NodeChildren(3).Title.Interpreter = 'none';
h.XLabel = [TTList{ti} ' scene selectivity index'];
colorbar off

customXLabels={};
for t=1:size(targ,2)
customXLabels{t} = jjnum2str(u(t),2); 
end
h.XDisplayLabels = customXLabels; 
h.CellLabelColor = 'none';
saveas(gca,[ROOT.Fig '\' field_name '_' TTList{ti} '.png'])
end
end