% display_variable

Initial_SWRFilter_common;
warning off
ROOT.Save = [ROOT.Mother '\Processed Data'];
ROOT.Rip0 = [ROOT.Mother '\Processed Data\ripples_mat\R0'];
ROOT.Rip = [ROOT.Mother '\Processed Data\ripples_mat\R3'];
ROOT.Rip4 = [ROOT.Mother '\Processed Data\ripples_mat\R4'];
ROOT.Fig3 = [ROOT.Mother '\Processed Data\ripples_mat\ProfilingSheet\R11_sub_field'];
ROOT.Units = [ROOT.Mother '\Processed Data\units_mat\U1'];
ROOT.Behav = [ROOT.Mother '\Processed Data\behavior_mat'];

dir = '';
ROOT.Fig = [ROOT.Fig3];
if ~exist(ROOT.Fig), mkdir(ROOT.Fig); end


Recording_region = readtable([ROOT.Info '\Recording_region_SWR.csv'],'ReadRowNames',true);
SessionList = readtable([ROOT.Info '\SessionList_SWR.xlsx'],'ReadRowNames',false);


%%
thisRegion = 'SUB';
thisRegion2 = 'SUB_field';
RipplesTable = readtable([ROOT.Save '\RipplesTable_' thisRegion2 '_forAnalysis_RDI.xlsx']);
% UnitsTable = readtable([ROOT.Units '\UnitsTable_filtered_' thisRegion2 '.xlsx']);
UnitsTable_A = readtable([ROOT.Units '\UnitsTable_' thisRegion '_forAnalysis.xlsx']);
UnitsTable_B = readtable([ROOT.Units '\UnitsTable_' thisRegion2 '_forAnalysis.xlsx']);

%%
nRDI_MF = RipplesTable.nRDI_MF;
nCells = RipplesTable.nTPPs;
x = round((1-nRDI_MF).*nCells); % 0~5까지의 랜덤한 인티저 수열 생성
% y = histcounts(x, 0:(max(x)+1)); % 각 숫자의 갯수 세기
% labels = string(y); % 갯수로 구성된 문자열 벡터 생성
p=pie(y, labels); 
title('# of single field cells')
set(gca,'FontSize',15,'FontWeight','bold')
t = findobj(p,'Type','text'); % 텍스트 객체 찾기
set(t,'FontSize',15,'FontWeight','bold'); % 텍스트 객체의 글꼴 크기와 굵기 변경

%% single field + univariate (pie chart)

nCells = RipplesTable.nTPPs;

figure;
sgtitle('# of single field + univariate cells')

subplot(2,3,1)
nRDI_UV = RipplesTable.nRDI_MF;
x = round((1-nRDI_UV).*nCells);
y = histcounts(x, 0:(max(x)+1)); % 각 숫자의 갯수 세기
labels = string(y); % 갯수로 구성된 문자열 벡터 생성
labels(strcmp(labels,'0'))="";
p=pie(y, labels); 
title('(single field only)')
caxis([0 11])

subplot(2,3,2)
nRDI_UV = RipplesTable.nRDI_hetero_L;
x = round((1-nRDI_UV).*nCells);
y = histcounts(x, 0:(max(x)+1)); % 각 숫자의 갯수 세기
labels = string(y); % 갯수로 구성된 문자열 벡터 생성
labels(strcmp(labels,'0'))="";
p=pie(y, labels); 
title('Left Scene')
caxis([0 11])

subplot(2,3,3)
nRDI_UV = RipplesTable.nRDI_hetero_R;
x = round((1-nRDI_UV).*nCells);
y = histcounts(x, 0:(max(x)+1)); % 각 숫자의 갯수 세기
labels = string(y); % 갯수로 구성된 문자열 벡터 생성
labels(strcmp(labels,'0'))="";
p=pie(y, labels); 
title('Right Scene')
caxis([0 11])

subplot(2,3,5)
nRDI_UV = RipplesTable.nRDI_hetero_C;
x = round((1-nRDI_UV).*nCells);
y = histcounts(x, 0:(max(x)+1)); % 각 숫자의 갯수 세기
labels = string(y); % 갯수로 구성된 문자열 벡터 생성
labels(strcmp(labels,'0'))="";
p=pie(y, labels); 
title('Choice')
caxis([0 11])

subplot(2,3,6)
nRDI_UV = RipplesTable.nRDI_hetero_SC;
x = round((1-nRDI_UV).*nCells);
y = histcounts(x, 0:(max(x)+1)); % 각 숫자의 갯수 세기
labels = string(y); % 갯수로 구성된 문자열 벡터 생성
labels(strcmp(labels,'0'))="";
p=pie(y, labels); 
title('Scene vs. Choice')
caxis([0 11])

colormap jet

%% single field + univariate (histogram)

nCells = RipplesTable.nTPPs;

figure;
sgtitle('# of single field + univariate cells')

subplot(2,3,1)
nRDI_UV = RipplesTable.nRDI_MF;
x = round((1-nRDI_UV).*nCells);
histogram(x, 'binedges',[0:12]); 
title('(single field only)')
xlabel('single field cells'); ylabel('ripples')

subplot(2,3,2)
nRDI_UV = RipplesTable.nRDI_hetero_L;
x = round((1-nRDI_UV).*nCells);
histogram(x, 'binedges',[0:12]); 
title('Left Scene')
xlabel('single+univariate cells'); ylabel('ripples')

subplot(2,3,3)
nRDI_UV = RipplesTable.nRDI_hetero_R;
x = round((1-nRDI_UV).*nCells);
histogram(x, 'binedges',[0:12]); 
title('Right Scene')
xlabel('single+univariate cells'); ylabel('ripples')

subplot(2,3,5)
nRDI_UV = RipplesTable.nRDI_hetero_C;
x = round((1-nRDI_UV).*nCells);
histogram(x, 'binedges',[0:12]); 
title('Choice')
xlabel('single+univariate cells'); ylabel('ripples')

subplot(2,3,6)
nRDI_UV = RipplesTable.nRDI_hetero_SC;
x = round((1-nRDI_UV).*nCells);
histogram(x, 'binedges',[0:12]); 
title('Scene vs. Choice')
xlabel('single+univariate cells'); ylabel('ripples')

%% get ratio of non-spatial ripple using univariate+SF

figure;
x = [RipplesTable.pRDI_L_UV, RipplesTable.pRDI_R_UV, RipplesTable.pRDI_C_UV];
x0 = x<0.05;

x1 = max(x0,[],2);

x2 = RipplesTable.DecodingP_all<0.05;

x3 = x1+2*x2;

id = min(isnan(x),[],2);

% x3(id)=[];

y = histcounts(x3, 0:(max(x3)+1)); % 각 숫자의 갯수 세기
y = y([1,3,4,2]);
labels = string(y); % 갯수로 구성된 문자열 벡터 생성
labels(strcmp(labels,'0'))="";
p=pie(y, labels); 



x = [RipplesTable.pRDI_L_UV_scuv, RipplesTable.pRDI_R_UV_scuv, RipplesTable.pRDI_C_UV_scuv];

%%
t = '_SC';
x = RipplesTable.(['nRDI_hetero' t]);
y = nanmean(UnitsTable_A.(['MultiVar' t]));

figure
subplot(1,2,1)
histogram((x))
line([y y], [0 120],'color','r')
ylabel('# of ripples')
title(['RDI_' t],'Interpreter','none')
subplot(1,2,2)
set(cdfplot((x)), 'LineWidth', 3); 
line([y y], [0 1],'color','r')
ylabel('ripple proportion')
title(['RDI_' t],'Interpreter','none')



%%
mRDI_SF = 1-RipplesTable.nRDI_MF;
nRDI_SF = round((1-nRDI_MF).*nCells);
pReplay = logical(RipplesTable.DecodingP_all<0.05);


pRDI = logical(RipplesTable.pRDI_L_SF<0.05 | RipplesTable.pRDI_R_SF<0.05 | RipplesTable.pRDI_C_SF<0.05);

figure

nbins = 10; % 구간의 개수
edges = linspace(min(nRDI_SF),max(nRDI_SF),nbins+1); % 구간의 경계값
counts = histcounts(nRDI_SF,edges); % 각 구간에 속하는 nRDI_SF의 개수
ratios = zeros(1,nbins); % 각 구간에서 pReplay가 true인 비율
for i = 1:nbins
    idx = (nRDI_SF >= edges(i)) & (nRDI_SF < edges(i+1)); % i번째 구간에 속하는 인덱스
    ratios(i) = mean(pRDI(idx)); % i번째 구간에서 pReplay가 true인 비율
end
bar(ratios) % 막대그래프 그리기
xticklabels(string(edges(1:end-1))) % x축 레이블 지정
xlabel('nRDI\_SF') % x축 이름 지정
ylabel('Ratio of pRDI') % y축 이름 지정
xlim([0 20])

%%
pRDI_R = (RipplesTable.pRDI_L_SF<0.05 | RipplesTable.pRDI_R_SF<0.05); % 첫 번째 조건에 해당하는 logical 배열
pRDI_C = (RipplesTable.pRDI_C_SF<0.05); % 두 번째 조건에 해당하는 logical 배열
pRDI = pRDI_R + pRDI_C*2; % 두 배열을 더하여 pRDI 구하기

nRDI_SF = categorical(nRDI_SF);
pRDI = categorical(pRDI);

pRDI2 = removecats(pRDI,'0'); % pRDI에서 0인 경우 제거
vals2 = categories(pRDI2); % pRDI2의 값
nvals2 = numel(vals2); % 값의 개수
counts2 = zeros(nvals2,ncats); % 각 값에 대해 각 카테고리에서 발생한 갯수
for i = 1:ncats
    idx = (nRDI_SF == cats(i)); % i번째 카테고리에 속하는 인덱스
    for j = 1:nvals2
        idy = (pRDI2 == vals2(j)); % j번째 값에 해당하는 인덱스
        counts2(j,i) = sum(idy(idx)); % j번째 값에 대해 i번째 카테고리에서 발생한 갯수
    end
end

total = sum(counts2); % 각 카테고리에서 pRDI가 0이 아닌 경우의 총 갯수
ratios2 = counts2 ./ total; % 총 갯수로 나누어 비율 구하기

b = bar(ratios([2 4 3],:)','stacked'); % 누적 막대그래프 그리기
xticklabels(cats) % x축 레이블 지정
xlabel('nRDI\_SF') % x축 이름 지정
ylabel('Ratio of pRDI') % y축 이름 지정
legend(vals2([1 3 2]),'Location','northoutside','orientation','horizontal') % 범례 표시
ylim([0 1])


xtips = b(1).XEndPoints; % 막대 끝 지점의 x좌표
ytips = sum(ratios2); % 막대 끝 지점의 y좌표
labels = string(counts2); % 갯수로 구성된 문자열 행렬 생성
labels = labels(:)'; % 문자열 행렬을 벡터로 변환
text(xtips,ytips,labels,'HorizontalAlignment','center',... % 막대 위에 텍스트 표시
    'VerticalAlignment','bottom')
