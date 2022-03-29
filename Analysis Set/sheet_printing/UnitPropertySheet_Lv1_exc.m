% ------------------ visual display
%
% units should be pixels. It's more handible.

sheet = figure; 
sheetPOS = [2 2 29.7 21.0];
% sheetPOS = [0 0 29.7 21.0];
set(gcf,'units','centimeters','position',sheetPOS,'Color', [1 1 1]);

szFONT = 10; 
txtFONT = 14;

%% clusterID & printed date
subplot(5,6,1:6);
set(gca,'units','normalized','position',[0.03 0.97 0.94 0.02]); axis off;
text(0.02, 0.5, ['Rat' thisRID '-S' thisSID '-TT' thisTTID '-cluster' thisCLID],'FontSize',txtFONT);
text(0.85, 0.5, ['printed on 25-Mar-2022'],'FontSize',szFONT);
% text(0.85, 0.5, ['printed on ' date],'FontSize',szFONT);

%% waveform average
aX = 10; aY = 590; aWidth = 125; aHeight = 125;

for ttRUN = 1:1:4
    subplot(5,6,ttRUN+6);
%     subplot(1,4,ttRUN);
    set(gca,'units','pixels','position',[aX aY aWidth aHeight]);
    aX = aX + 135;
    
    % wait bar
    w = waitbar(0,['Please wait... (' num2str(ttRUN) '/4)']);
    for clRUN = 1:1:size(thisEpochCLAPrand, 3)
        plot(gca,thisEpochCLAPrand(:,ttRUN,clRUN), 'Color', [.8 .8 .8]);
        waitbar(clRUN/size(thisEpochCLAPrand, 3));
        hold on;
    end
    plot(gca,thisMEANAP(:,ttRUN), 'Color', [0 0 0]);
    axis([0 32 min(min(min(thisEpochCLAPrand)))-100 max(max(max(thisEpochCLAPrand)))+100]); axis off;
    title([jjnum2str(max(thisMEANAP(:,ttRUN)),1) ' \muv'],'FontSize',8); %ylabel(['max ' num2str(max(max(max(thisEpochCLAP)))) 'mV']);
    hold off;
    close(w);
end


%% time vs. max height
subplot(5,6,11);
set(gca,'units','pixels','position',[585 600 270 130]);

height_range = max(squeeze(thisEpochCLAP(:,max_channel,:)))-min(squeeze(thisEpochCLAP(:,max_channel,:)));
plot(thisEpochCLTS/1000000, height_range, 'p','MarkerEdgeColor','k', 'MarkerSize', 1);
axis([thisEpochCLTS(1)/1000000, thisEpochCLTS(end)/1000000, 0, max(height_range)+(max(height_range)-min(height_range))*0.15]);
title('Time vs Max peak','FontWeight','Bold','FontSize',szFONT);
xlabel(['Time (s)']);

%% log ISI
subplot(5,6,11);
set(gca,'units','pixels','position',[900 580 145 150]);
hold on;
isiHIST(400) = max(isiHIST)-1;
bar(isiHIST);	%plot(ones(1, max(isiHIST) + 1) .* find(histEDGE == 0), 0:1:max(isiHIST), 'r:', 'LineWidth', szLINE);
text(min(find(isiHIST == min(max(isiHIST)))), min(max(isiHIST)), [jjnum2str((10^histEDGE(min(find(isiHIST == min(max(isiHIST)))))) / 1000, 1) ' ms']); LogISIPEAKTIME = (10^histEDGE(min(find(isiHIST == min(max(isiHIST)))))) / 1000;
% text(400,max(isiHIST),'1ms');
hold off;
title(['log ISI'],'FontWeight','Bold','FontSize',szFONT); xlabel(['time (ms)']); axis tight; 
set(gca, 'FontSize', szFONT, 'XLim', [350 size(histEDGE, 2)], 'XTick', 350:((size(histEDGE, 2) - 350) / 2):size(histEDGE, 2), 'XTickLabel', {['.31'], ['55'], ['10000']});
% title(['log ISI']); xlabel(['time (ms)']); axis tight; set(gca, 'FontSize', szFONT, 'XLim', [0 size(histEDGE, 2)], 'XTick', 0:((size(histEDGE, 2)) / 2):size(histEDGE, 2), 'XTickLabel', {['.31'], ['55'], ['10000']});

%% autocorrelation
subplot(5,6,29);
set(gca,'units','pixels','position',[900 310 180 180]);
bar(corrXlabel, correlogram);
title(['AutoCorrelogram'],'FontWeight','Bold','FontSize',szFONT); set(gca, 'XLim', [min(corrXlabel) max(corrXlabel)], 'YLim', [0 max(log10(correlogram))]); xlabel(['Time (ms)']); axis tight;

%% skaggs firing rate map
% firing map bin size = (45x30) x 7.4
imCOL = 222; imROW = 333;


subplot(5,6,13);

set(gca,'units','pixels','position',[15 240  imCOL imROW]);  %axis off;
set(gca, 'YDir', 'rev', 'nextplot', 'add','XLim', [0 300 / thisFRMapSCALE], 'YLim', [0 450 / thisFRMapSCALE]);
if strcmp(exper,'JS')
    skaggsMap_resized = flip(skaggsMap(:,1));
    set(gca,'units','pixels','position',[15 240  imCOL imROW]);  %axis off;
    set(gca, 'nextplot', 'add','XLim', [-2 4], 'YLim', [0 size(skaggsMap_resized,1)]);
else
skaggsMap_resized = skaggsMap(6:50,21:50);
end
[r c] = size(skaggsMap_resized);
max_ratemap = max(max(skaggsMap_resized));
skaggsMap_resized = max_ratemap - skaggsMap_resized;

thisAlphaZ = skaggsMap_resized; thisAlphaZ(isnan(skaggsMap_resized)) = 0; thisAlphaZ(~isnan(skaggsMap_resized)) = 1;
colormap gray;
skaggsMap_alpha = skaggsMap_resized;
skaggsMap_alpha(isnan(skaggsMap_alpha)) = inf; 
imagesc(skaggsMap_alpha); % alpha(thisAlphaZ); 
% title(['Skaggs firing rate map'],'FontWeight','Bold','FontSize',szFONT);
set(gca,'xtick',[],'ytick',[]); axis off;

% rate map scale
subplot(4,6,23);
set(gca,'units','pixels','position',[50 190 170 20]);
set(gca, 'XLim', [0 10], 'YLim', [0 10]); axis off;

ratemap_scale = [];
ratemap_scale(:,1) = [max_ratemap : -max_ratemap/50 : 0];
ratemap_scale = transpose(ratemap_scale);

imagesc(ratemap_scale);
set(gca,'xtick',[1 50],'XTickLabel',{['0'] [jjnum2str(onmazeMaxFR,2) ' Hz']},'ytick',[]);


%% raw map
szDOT = 2;
colSPK = [0 0 0]; % changed from RGB to gray scale

xlim = [200 500];
ylim = [50 500];
if strcmp(exper,'JS')
    xlim = [-2 4];
    ylim = [0 1000];
end
    
imCOL = 222; % imCOL = x range * 0.74.
imROW = 300; % imCOL = x range * 0.74.

% occJPG is plotted with dots of [0.7 0.7 0.7] color.

subplot(5,6,15);
hold on;
set(gca,'units','pixels','position',[300 240 imCOL imROW]); 
scatter(Pos.x(in) - xlim(1), Pos.y(in)  - ylim(1)+1,10,[.8 .8 .8],'filled')
plot(Spk.x_spk(inspk) - xlim(1), Spk.y_spk(inspk) - ylim(1)+1, '.', 'MarkerSize', szDOT, 'color', colSPK);

title('raw spk map','FontWeight','Bold','FontSize',szFONT)

if strcmp(exper,'JS')
set(gca, 'XLim', [0 xlim(2)-xlim(1)], 'YLim', [0 ylim(2)-ylim(1)],'xtick',[],'ytick',[]);
else
   set(gca, 'YDir', 'rev','XLim', [0 xlim(2)-xlim(1)], 'YLim', [0 ylim(2)-ylim(1)],'xtick',[],'ytick',[]); 
end
axis off;


if max(Pos.trial_context)==4
    subplot('position',[.52 .51 .1 .15]);
    DrawRawSpkMap(Pos,Spk,Title,1,xlim,ylim,in,inspk,exper)
    subplot('position',[.64 .51 .1 .15]);
    DrawRawSpkMap(Pos,Spk,Title,2,xlim,ylim,in,inspk,exper)
    subplot('position',[.52 .32 .1 .15]);
    DrawRawSpkMap(Pos,Spk,Title,3,xlim,ylim,in,inspk,exper) 
    subplot('position',[.64 .32 .1 .15]);
    DrawRawSpkMap(Pos,Spk,Title,4,xlim,ylim,in,inspk,exper)
end

if max(Pos.trial_context)==2
    subplot('position',[.52 .32 .1 .32]);
    DrawRawSpkMap(Pos,Spk,Title,1,xlim,ylim,in,inspk,exper)
    subplot('position',[.64 .32 .1 .32]);
    DrawRawSpkMap(Pos,Spk,Title,2,xlim,ylim,in,inspk,exper)
end
%% text
txtX1 = 1; txtX2 = 6; txtValueX1 = 4; txtValueX2 = 9;
txtINIY = 9; txtADJY = 1.6;
txtFONT = 12; 

subplot(4,6,24);
set(gca,'units','pixels','position',[240 10 840 200]);
set(gca, 'XLim', [0 10], 'YLim', [0 10]); axis off;

text(txtX1, txtINIY - txtADJY * 0, ['Spike height (from baseline)'], 'FontSize', txtFONT);
text(txtX1, txtINIY - txtADJY * 1, ['Spike height (peak-to-valley)'], 'FontSize', txtFONT);
text(txtX1, txtINIY - txtADJY * 2, ['Spike width (peak-to-valley)'], 'FontSize', txtFONT);
text(txtX1, txtINIY - txtADJY * 3, ['Session Avg. FR'], 'FontSize', txtFONT);
text(txtX1, txtINIY - txtADJY * 4, ['On-map Avg. FR'], 'FontSize', txtFONT);
text(txtX1, txtINIY - txtADJY * 5, ['On-map Max FR'], 'FontSize', txtFONT);

text(txtValueX1, txtINIY - txtADJY * 0, [jjnum2str(max_peak,1) ' (\muv)'], 'FontSize', txtFONT);
text(txtValueX1, txtINIY - txtADJY * 1, [jjnum2str(max_amp,1) ' (\muv)'], 'FontSize', txtFONT);
text(txtValueX1, txtINIY - txtADJY * 2, [jjnum2str(max_width,1) ' (\mus)'], 'FontSize', txtFONT);
text(txtValueX1, txtINIY - txtADJY * 3, [jjnum2str(FRRate,2) ' (Hz)'], 'FontSize', txtFONT);
text(txtValueX1, txtINIY - txtADJY * 4, [jjnum2str(onmazeAvgFR,2) ' (Hz)'], 'FontSize', txtFONT);
text(txtValueX1, txtINIY - txtADJY * 5, [jjnum2str(onmazeMaxFR,2) ' (Hz)'], 'FontSize', txtFONT);

text(txtX2, txtINIY - txtADJY * 0, ['Spatial Information Score' ], 'FontSize', txtFONT);
text(txtX2, txtINIY - txtADJY * 1, ['# of spks total'], 'FontSize', txtFONT);
text(txtX2, txtINIY - txtADJY * 2, ['# of spks in boundary' ], 'FontSize', txtFONT);
text(txtX2, txtINIY - txtADJY * 3, ['Spikes within refractory period' ], 'FontSize', txtFONT);
% text(txtX2, txtINIY - txtADJY * 4, ['L-Ratio'], 'FontSize', txtFONT);
% text(txtX2, txtINIY - txtADJY * 5, ['Isolation Distance'], 'FontSize', txtFONT);

text(txtValueX2, txtINIY - txtADJY * 0, [jjnum2str(SpaInfoScore,2) ' (bit/spk)'], 'FontSize', txtFONT);
text(txtValueX2, txtINIY - txtADJY * 1, [num2str(nSPKS)], 'FontSize', txtFONT);
text(txtValueX2, txtINIY - txtADJY * 2, [num2str(nSPKS_in) ], 'FontSize', txtFONT);
text(txtValueX2, txtINIY - txtADJY * 3, [jjnum2str(withinREFRACPortion,2) ' %'], 'FontSize', txtFONT);
% text(txtValueX2, txtINIY - txtADJY * 4, [jjnum2str(LRATIO,2) ' '], 'FontSize', txtFONT);
% text(txtValueX2, txtINIY - txtADJY * 5, [jjnum2str(ISODIST,2) ' '], 'FontSize', txtFONT);
%%
cd(ROOT.Save)
saveImage(sheet, ['rat' clusterID '.png'], 'centimeters', sheetPOS);

%%
function DrawRawSpkMap(Pos,Spk,Title,num,xlim,ylim,in,inspk,exper)
if ~isempty(Title{num})

scatter(Pos.x(in) - xlim(1), Pos.y(in)  - ylim(1)+1,10,[.8 .8 .8],'filled')
hold on
plot(Spk.x_spk(Spk.cont_spk(:,num)) - xlim(1), Spk.y_spk(Spk.cont_spk(:,num)) - ylim(1)+1, '.', 'MarkerSize', 2, 'color', 'k');
title(Title(num),'FontWeight','Bold','FontSize',10)

if strcmp(exper,'JS')
set(gca,'XLim', [0 xlim(2)-xlim(1)], 'YLim', [0 ylim(2)-ylim(1)],'xtick',[],'ytick',[]);
else
   set(gca, 'YDir', 'rev','XLim', [0 xlim(2)-xlim(1)], 'YLim', [0 ylim(2)-ylim(1)],'xtick',[],'ytick',[]); 
end
end
axis off;

end
