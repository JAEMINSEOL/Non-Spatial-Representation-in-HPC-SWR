Initial_SWRFilter_common;
warning off


FRbin=[]; thisUnit=[]; UnitSum_Ripples=table;
thisRegion = 'CA1';
crit_cells=.5;
ReactProb = [];
UnitsTable_all=[];

%%

UnitsTable_all=[];
RipplesTable_all = [];
%%
for clRUN=1:5
    %%
    thisSID.full = ['561-0' num2str(clRUN)];
    load([ROOT.Save '\' thisSID.full '_' thisRegion '.mat'])
    Behav = load([ROOT.Raw.Mother '\rat' thisSID.full(1:3) '\rat' thisSID.full '\' 'ParsedPosition.mat']);
    [Behav, trial_vector] = MakeTrialVector(Behav);
    %%
    thisUnit=[];
    for clRip = 1:size(RipplesTable,1)
        t = floor(RipplesTable.Trial(clRip));
        if t==0, Cxt=0;
        else
            Cxt = Behav.trial_context(t);
        end
        
        if RipplesTable.Trad(clRip)==0
            continue;
        end
        
        
        RipID = cell2mat(RipplesTable.RippleID(clRip));
        thisRippleStruct = RipplesStruct.(thisRegion).(['ripple' RipID(end-3:end)]);
        thisUnit_temp = thisRippleStruct.SpkList.first(:,1:8);
        thisUnit_temp(:,1) = thisUnit_temp(:,1) ./ RipplesTable.RippleDuration(clRip);
        thisUnit_temp = [thisUnit_temp(:,1), repmat([RipplesTable.Trad(clRip)],size(thisUnit_temp,1),1), repmat([561 clRUN],size(thisUnit_temp,1),1),thisUnit_temp(:,2:end)...
            ];
        
        thisUnit_temp = [thisUnit_temp, repmat(double(Cxt),size(thisUnit_temp,1),1),repmat(clRip,size(thisUnit_temp,1),1)];
        
        id = thisUnit_temp(:,6)>=crit_cells;
        thisUnit = [thisUnit;thisUnit_temp];
    end
    %%
    for clUnit = 1:size(UnitsTable,1)
        clusterID = [thisSID.full '-' num2str(UnitsTable.TT(clUnit)) '-' jmnum2str(UnitsTable.Unit(clUnit),2)];
        if UnitsTable.Type(clUnit)<=0, continue; end
        if exist([ROOT.Raw.Map '\rat' clusterID '.mat'])
            Spk_processed = load([ROOT.Raw.Map '\rat' clusterID '.mat']);
            Dv = Spk_processed.stem_end_index;
            FRMap = [];
            for i=1:5
                FRMap(i,:) = Spk_processed.skaggsMap1D{i};
            end
            
            FRMap(6,:) = Spk_processed.skaggsMap_left1D;  
            FRMap(7,:) = Spk_processed.skaggsMap_right1D; 
            
            %      if Spk_processed.onmazeAvgFR1D(1)<1
            %          continue;
            %      end
            
            
            
            %%
            if isnan(UnitsTable.RDI_LR(clUnit)), continue; end
            
            % reactivation probability
            cxt=[];
            for i=1:size(RipplesTable.Trial,1)
                if RipplesTable.Trial(i)<1
                    cxt(i,:)=0;
                else
                    cxt(i,:) = Behav.trial_context(floor(RipplesTable.Trial(i)));
                end
            end
            nRip=[];
            
            nRip(1,1) = sum(RipplesTable.Filter>0);
            nRip(2,1) = sum(RipplesTable.Filter>0 & cxt==1);
            nRip(3,1) = sum(RipplesTable.Filter>0 & cxt==3);
            nRip(4,1) = sum(RipplesTable.Filter>0 & cxt==2);
            nRip(5,1) = sum(RipplesTable.Filter>0 & cxt==4);
            
            nRip(6,1) = sum(RipplesTable.Filter>0 & (cxt==1|cxt==3));
            nRip(7,1) = sum(RipplesTable.Filter>0 & (cxt==2|cxt==4));
            
            
            
            nSpk =[sum(cell2mat(UnitsTable.NumRipples(clUnit,:))); cell2mat(UnitsTable.NumRipples(clUnit,[1,3,2,4]))';...
                sum(cell2mat(UnitsTable.NumRipples(clUnit,1:2:3)));sum(cell2mat(UnitsTable.NumRipples(clUnit,2:2:4)))];
            
            if sum(nSpk(1,:))==0, close all; continue; end
            
            v = (nSpk(:,1)./nRip(:,1));
            
            
            %%
            FigPos.Background = [100 100 1400 800];
            FigPos.LinFRMaps = [.05 .65 .2 .2];
            FigPos.TextBox = [.05 .9 .4 .1];
            FigPos.ReactProp = [.06 .05 .11 .2];
            FigPos.RDI_bar = [.19 .32 .03 .2];
            FigPos.NormReact = [.54 .56 .18 .04];
            FigPos.ReactProbScatter = [.61 .48 .14 .14];
            
            clist = {'#E72416','#4FBFD6','#711419','#4F63AD','#FF0000','#0000FF'};
            
            figure('position',FigPos.Background,'color','w')

            annotation('textbox',FigPos.TextBox,'String',sprintf(['rat' clusterID '\n' ...
                '                RDI:   ZB=' jjnum2str(UnitsTable.RDI_ZB(clUnit),2) ',            PM= ' jjnum2str(UnitsTable.RDI_PM(clUnit),2) ...
                ',         LR= ' jjnum2str(UnitsTable.RDI_LR(clUnit),2)...
                '\nRip. Part. Rate: Z=' jjnum2str(v(1),2) ', B=' jjnum2str(v(2),2) ', P=' jjnum2str(v(3),2) ', M=' jjnum2str(v(4),2) ... 
                ', L=' jjnum2str(v(5),2) ', R=' jjnum2str(v(6),2) ... 
                '\nSI = ' jjnum2str(Spk_processed.SpaInfoScore1D(1),2)]),...
                'EdgeColor','none','FontSize',10)
            
            
            subplot('position',FigPos.LinFRMaps)
            hold on
            plot(FRMap(1,:),'color','k','linewidth',2)
            s1 = ScatMaxFR(FRMap(1,:));
            SetBasicPropFR(FRMap,Dv)
            title('overall','fontweight','b')
            
            subplot('position',FigPos.LinFRMaps+[.23 0 0 .05])
            hold on
            l1 = plot(FRMap(2,:),'color',clist{1},'linewidth',1);
            s1 = ScatMaxFR(FRMap(2,:));
            l2 = plot(FRMap(3,:),'color',clist{3},'linewidth',1);
            s2 = ScatMaxFR(FRMap(3,:));
            SetBasicPropFR(FRMap,Dv)
            title('left choice','fontweight','b')
            legend([l1 l2],{'Zebra','Bamboo'},'location','northoutside','orientation','horizontal')
            
            subplot('position',FigPos.LinFRMaps+[.46 0 0 .05])
            hold on
            l1 = plot(FRMap(4,:),'color',clist{2},'linewidth',1);
            s1 = ScatMaxFR(FRMap(4,:));
            l2 = plot(FRMap(5,:),'color',clist{4},'linewidth',1);
            s2 = ScatMaxFR(FRMap(5,:));
            SetBasicPropFR(FRMap,Dv)
            title('right choice','fontweight','b')
            legend([l1 l2],{'Pebbles','Mountain'},'location','northoutside','orientation','horizontal')
            
            subplot('position',FigPos.LinFRMaps+[.69 0 0 .05])
            hold on
            l1 = plot(FRMap(6,:),'color',hex2rgb(clist(5)),'linewidth',1);
            s1 = ScatMaxFR(FRMap(6,:));
            l2 = plot(FRMap(7,:),'color',hex2rgb(clist(6)),'linewidth',1);
            s2 = ScatMaxFR(FRMap(7,:));
            SetBasicPropFR(FRMap,Dv)
            title('left vs. right choice','fontweight','b')
            legend([l1 l2],{'Left','Right'},'location','northoutside','orientation','horizontal')
            

             m = max(FRMap,[],2);
             n = nanmean(FRMap,2);
             
             subplot('position',FigPos.ReactProp+[0 .3 -.03 0])
            b= DrawRprPlot_2items(m(1),'k','w','Max. Firing Rate (Hz)',max(m)*1.3);
             xticks([1]); xticklabels({'All'});
             set(b,'facealpha',.5)
%               subplot('position',FigPos.ReactProp+[.1 .3 -.03 0])
             b=DrawRprPlot_2items(n(1),'k','w','Mean Firing Rate (Hz)',max(m)*1.3);
             xticks([1]); xticklabels({'All'});
             
             subplot('position',FigPos.ReactProp+[.23 .3 -.03 0])
             b=DrawRprPlot_2items(m(2:3),hex2rgb(clist(1)),hex2rgb(clist(3)),'Max. Firing Rate (Hz)',max(m)*1.3);
             xticks([1:2]); xticklabels({'Z','B'})
             set(b,'facealpha',.5)
%              subplot('position',FigPos.ReactProp+[.34 .3 -.03 0])
             b=DrawRprPlot_2items(n(2:3),hex2rgb(clist(1)),hex2rgb(clist(3)),'Max/Mean Firing Rate (Hz)',max(m)*1.3);
             xticks([1:2]); xticklabels({'Z','B'})
             
             subplot('position',FigPos.ReactProp+[.46 .3 -.03 0])
             b=DrawRprPlot_2items(m(4:5),hex2rgb(clist(2)),hex2rgb(clist(4)),'Max. Firing Rate (Hz)',max(m)*1.3);
             xticks([1:2]); xticklabels({'P','M'})
             set(b,'facealpha',.5)
%              subplot('position',FigPos.ReactProp+[.57 .3 -.03 0])
             b=DrawRprPlot_2items(n(4:5),hex2rgb(clist(2)),hex2rgb(clist(4)),'Max/Mean Firing Rate (Hz)',max(m)*1.3);
             xticks([1:2]); xticklabels({'P','M'})
             
             subplot('position',FigPos.ReactProp+[.69 .3 -.03 0])
             b=DrawRprPlot_2items(m(6:7),hex2rgb(clist(5)),hex2rgb(clist(6)),'Max. Firing Rate (Hz)',max(m)*1.3);
             xticks([1:2]); xticklabels({'L','R'})
             set(b,'facealpha',.5)
%              subplot('position',FigPos.ReactProp+[.8 .3 -.03 0])
             b=DrawRprPlot_2items(n(6:7),hex2rgb(clist(5)),hex2rgb(clist(6)),'Max/Mean Firing Rate (Hz)',max(m)*1.3);
             xticks([1:2]); xticklabels({'L','R'})
             
             subplot('position',FigPos.ReactProp)
             hold on
             %             b1 = barh(nRip(1,:));
             %             b1(1).FaceColor = 'w';
             %             b2 = barh(nSpk(1,:));
             %            b2(1).FaceColor = 'k';
             
             v = (nSpk(:,1)./nRip(:,1));
             b=DrawRprPlot_2items(v([1],:),'k','w','Ripple Participation Rate',.6);
             xticks([1]); xticklabels({'All'});
%             title('Ripple Participation Rate')
            
            subplot('position',FigPos.ReactProp+[.25 0 0 0])
           b=DrawRprPlot_2items(v([2,3],:),hex2rgb(clist(1)),hex2rgb(clist(3)),'Ripple Participation Rate',max(v)*1.3);
            xticks([1:2]); xticklabels({'Z','B'})
            
            subplot('position',FigPos.ReactProp+[.48 0 0 0])
            b=DrawRprPlot_2items(v([4,5],:),hex2rgb(clist(2)),hex2rgb(clist(4)),'Ripple Participation Rate',max(v)*1.3);
            xticks([1:2]); xticklabels({'P','M'})
            
            subplot('position',FigPos.ReactProp+[.71 0 0 0])
            DrawRprPlot_2items(v([6,7],:),hex2rgb(clist(5)),hex2rgb(clist(6)),'Ripple Participation Rate',max(v)*1.3);
            xticks([1:2]); xticklabels({'L','R'})
            
            
            subplot('position', FigPos.RDI_bar+[.18 .08 .08 -.15])
            r = UnitsTable.RDI_ZB(clUnit);
            if r>0, c=clist{1}; elseif r<0, c=clist{3}; else, c='#000000'; end
            m=max(1.2,abs(r));
            DrawHLine(-m,m,.8,'vertical')
            q=quiver(0,0,r,0);
            q.Color = hex2rgb(c); q.LineWidth=2;
            scatter(r,0,40,hex2rgb(c),'filled')
            text(1,-1,sprintf(['Z']),'fontweight','b')
            text(-1,-1,sprintf(['P']),'fontweight','b')
            title(['RDI_{SCN-L}=' jjnum2str(r,2)])
            
            subplot('position', FigPos.RDI_bar+[.41 .08 .08 -.15])
            r = UnitsTable.RDI_PM(clUnit);
            if r>0, c=clist{2}; elseif r<0, c=clist{4}; else, c='#000000'; end
            m=max(1.2,abs(r));
            DrawHLine(-m,m,.8,'vertical')
            q=quiver(0,0,r,0);
            q.Color = hex2rgb(c); q.LineWidth=2;
            scatter(r,0,40,hex2rgb(c),'filled')
            text(1,-1,sprintf(['P']),'fontweight','b')
            text(-1,-1,sprintf(['M']),'fontweight','b')
            title(['RDI_{SCN-R}=' jjnum2str(r,2)])
            
            subplot('position', FigPos.RDI_bar+[.64 .08 .08 -.15])
            r = UnitsTable.RDI_LR(clUnit);
            if r>0, c=clist{5}; elseif r<0, c=clist{6}; else, c='#000000'; end
            m=max(1.2,abs(r));
            DrawHLine(-m,m,.8,'vertical')
            q=quiver(0,0,r,0);
            q.Color = hex2rgb(c); q.LineWidth=2;
            scatter(r,0,40,hex2rgb(c),'filled')
            text(1,-1,sprintf(['L']),'fontweight','b')
            text(-1,-1,sprintf(['R']),'fontweight','b')
            title(['RDI_C=' jjnum2str(r,2)])

            
            annotation('line',[.255 .255],[.05 .9],'linestyle','--')
            annotation('line',[.485 .485],[.05 .9],'linestyle','--')
            annotation('line',[.715 .715],[.05 .9],'linestyle','--')
            %% Save image
            num='13U';
            cd([ROOT.Save '\Plot\IndivUnits\'])
            mkdir(num)
            saveas(gca,[ROOT.Save '\Plot\IndivUnits\' num '\' num '-rat' clusterID '.png'])
            disp([clusterID ' is saved!'])
            close all
        end
        
        
    end

end

function s = ScatMaxFR(FRMap)
[n,m] = max(FRMap);
s=scatter(m,n,20,'r','filled');
text(m*1.1,n,[jjnum2str(n,2) 'Hz'])
end

function SetBasicPropFR(FRMap,Dv)
line([Dv Dv],[0 max(max(FRMap))],'color','r','linewidth',1.5)
set(gca,'ylim',[0 max(max(FRMap))*1.1])
xticks([1 Dv 45]); xticklabels({'Startbox','Dv','Fdwell'});
ylabel('FR (Hz)')
end

function DrawReactRaster(v,i)


if i==0
    % id = v(:,2)==0.5;
    %      r=rasterplot_JM(ones(sum(id),1),v(id,1),[0 0 0],.7);
    id = logical(ones(size(v,1),1));
    r=rasterplot_JM(ones(sum(id),1),v(id,1),[0 0 0],.7);
    %      id = v(:,2)==1;
    %      r=rasterplot_JM(ones(sum(id),1),v(id,1),[1 0 0],.7);
    
else
    id = v(:,3)==i;
    %     id = v(:,2)==0.5 & v(:,3)==i;
    r=rasterplot_JM(ones(sum(id),1),v(id,1),[0 0 0],.7);
    %      id = v(:,2)==1& v(:,3)==i;
    %      r=rasterplot_JM(ones(sum(id),1),v(id,1),[1 0 0],.7);
end

axis off
line([0 1],[.5 .5],'color','k')
line([0 0],[.25 .75],'color','k')
line([1 1],[.25 .75],'color','k')
text(0,.2,'0')
text(1,.2,'1')

end

function ReactScatter(v2,xl,yl,ti)
scatter(v2(1,1),v2(2,1),20,'k','*')
% scatter(v2(1,2),v2(2,2),20,'r','filled')
% scatter(v2(1,3),v2(2,3),20,'k','filled')
xlim([0 1])
ylim([0 1])
xlabel(xl)
ylabel(yl)
title(ti,'fontweight','b')
line([0 1], [.5 .5], 'color','k','linestyle','--')
line( [.5 .5], [0 1],'color','k','linestyle','--')
% legend({'all','seq','nseq'},'location','eastoutside')
end

function DrawHLine(Llim, Ulim, h,rotation)
hold on
if strcmp(rotation,'horizontal')
    line([0 0],[Llim Ulim],  'color','k')
    line([h*(-1) h],[Llim Llim],  'color','k')
    line([h*(-1) h],[Ulim Ulim],  'color','k')
    line([h*(-.7) h*(.7)],[0 0],  'color','k')
    text(-h-1,Ulim,num2str(Ulim))
    text(-h-1,Llim,num2str(Llim))
    ylim([-2 2.5])
else
    line([Llim Ulim], [0 0], 'color','k')
    line([Llim Llim], [h*(-1) h], 'color','k')
    line([Ulim Ulim], [h*(-1) h], 'color','k')
    line([0 0],[h*(-.7) h*(.7)],  'color','k')
    xlim([-2 2.5])
        text(Ulim,-h-1,num2str(Ulim))
    text(Llim,-h-1,num2str(Llim))
    set(gca, 'XDir','reverse')
end

axis off
end

function b=DrawRprPlot_2items(v,c1,c2,tit,m)
hold on
%             b1 = barh(nRip([i,j],:));
%             b1(1).FaceColor = 'w';
if length(v)==2
            b = bar(v);
            b.FaceColor = 'flat';
            b.CData = [c1; c2];
            ytips1 = b(1).YEndPoints*1.1;
            xtips1 = b(1).XEndPoints;
            text(xtips1,ytips1,{jjnum2str(v(1),2),jjnum2str(v(2),2)},'HorizontalAlignment','center','VerticalAlignment','bottom')
            
else
              b = bar(v(1));
            b(1).FaceColor = c1;
            ytips1 = b(1).YEndPoints*1.1;
            xtips1 = b(1).XEndPoints;
            text(xtips1,ytips1,{jjnum2str(v(1),2)},'HorizontalAlignment','center','VerticalAlignment','bottom')
            
end
            
             ylabel(tit)
             title(tit)
%             axis ij
            set(gca,'fontweight','b')
            ylim([0 m])
end
