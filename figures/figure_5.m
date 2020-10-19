%
% figure_5
% Part of the JoNmusic2020 code (github.com/octaveEtard/JoNmusic2020)
% Author: Octave Etard, 2020
%
% Plot Figure 5.
%
% 1 - Load behavioural results
% 2 - Run statistical comparisons
% 2 - Plot.
%
%% Load results
maxRT = 2; % in s
% where are the results located
saveFolder = JoNmusic2020.getPath('behav');
saveName = sprintf('clickPerformance_RT_%.1f.mat', maxRT);


d = load(fullfile(saveFolder, saveName));
instrumentOrder = d.instrumentOrder;
% True/False Positve Rate
TPR = d.TPR;
FPR = d.FPR;
nSub = numel(d.SID);

%%
% Compute d-prime
dp = norminv(TPR)-norminv(FPR);
% Wilcoxon paired signed rank test
pval_dp = signrank(dp(:,1),dp(:,2)); % d-prime when attending guitar vs piano

% Comparing T and F positive rates between the 2 instruments
pval_TF = zeros(2,1);
% Wilcoxon paired signed rank test
pval_TF(1) = signrank(TPR(:,1),TPR(:,2)); % T guitar vs T piano
pval_TF(2) = signrank(FPR(:,1),FPR(:,2)); % F guitar vs F piano
% Comparing T to F positive rates for each instrument
pval_GP = zeros(2,1);
% Wilcoxon paired signed rank test
pval_GP(1) = signrank(TPR(:,1),FPR(:,1)); % T vs F for guitar
pval_GP(2) = signrank(TPR(:,2),FPR(:,2)); % T vs F for piano

% FDR correction
pval_FDR = fdr([pval_TF;pval_GP]);
pval_TF = pval_FDR(1:2);
pval_GP = pval_FDR(3:4);

% significance level
alpha0 = 5/100;

% print some values in terminal
m = mean(dp,1);
fprintf('Average d'': %.1f %s, %.1f %s ; pval = %.1e\n\n',m(1),instrumentOrder{1},...
    m(2),instrumentOrder{2},pval_dp);
fprintf('TPR, %s vs %s, pval = %.1e\n',instrumentOrder{:},pval_TF(1));
fprintf('FPR, %s vs %s, pval = %.1e\n',instrumentOrder{:},pval_TF(2));
fprintf('%s, TPR vs FPR, pval = %.1e\n',instrumentOrder{1},pval_GP(1));
fprintf('%s, TPR vs FPR, pval = %.1e\n',instrumentOrder{2},pval_GP(2));


%% Plot
fts = 11;
lwd = 1;
mks = 5.5;

% capitalise
instrumentOrder = cellfun(@(s) sprintf('%s%s',upper(s(1)),s(2:end)),...
    instrumentOrder,'UniformOutput',false);
labels = {'True positive rate','False positive rate'};
topLabels = {'A','B','C','D'};

% jitter scatter plots in the x-direction
maxJitter = 0.15;

col = JoNmusic2020.plotStyleArgs(); % colors

lwd_sub = 0.6;
markerLineSpacing = 0.1;

nRow = 2;
nCol = 2;

nPlot = nRow * nCol;
jitter = linspace(-maxJitter,maxJitter,nSub);
% jitter = 2*maxJitter*rand(1,nSub)-maxJitter;
axx = gobjects(nPlot,1);

fig = figure;

% --- ROC curve
scatterObj = gobjects(2,1); % for easy labelling

axx(1) = subplot(nRow,nCol,1); hold on;
ax = axx(1);

plot([0,1],[0,1],'Color',col(3,:),'LineWidth',lwd_sub); % 1-1 line
for iInstru = 1:2
    scatterObj(iInstru) = plot(FPR(:,iInstru),TPR(:,iInstru),'o','MarkerSize',mks,'Color',col(iInstru,:));
end

ax.XAxis.Limits = [0,1];
ax.YAxis.Limits = [0,1];

ax.DataAspectRatio = [1,1,1];

legend(ax, scatterObj, instrumentOrder, 'box', 'off',...
    'Location', 'southeast',...
    'FontSize', fts);

ax.XAxis.Label.String = labels{2};
ax.YAxis.Label.String = labels{1};

% --- d-prime plot
yLimits = [0,4.5];

% significance text parameters
y0 = 4.2; % height (data scale) 
dy = yLimits(end)*2.5e-2; % vertical bits length
compLineArgs = {'Color','k','Linewidth',1}; % line parameters

axx(2) = subplot(nRow,nCol,2); hold on;
ax = axx(2);

% jitter points a little along the x axis
plot([1+markerLineSpacing; 2-markerLineSpacing] + jitter, dp', 'Color', col(3,:),'LineWidth',lwd_sub);
plot(1 + jitter, dp(:,1), 'o','MarkerSize',mks,'Color',col(1,:));
plot(2 + jitter, dp(:,2), 'o','MarkerSize',mks,'Color',col(2,:));

ax.XAxis.Limits = [1,2] + 2*[-maxJitter,maxJitter];
ax.XAxis.TickValues = [1, 2];
ax.XAxis.TickLabels = instrumentOrder;
ax.YAxis.Limits = yLimits;
ax.YAxis.Label.String = 'd-prime';
%
if pval_dp <= alpha0
    pltools.addComparisonLine(ax,1:2,y0,dy,compLineArgs{:});
    txt = pltools.equalText('p',pval_dp,1);
    text(1.5,y0 + dy,txt,...
        'Units','data','HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontSize',round(0.8*fts))
end


% --- TPR and FPR plots
yLimits = [1,0.6];
% significance text parameters
y0 = [1,0.55];
dy = yLimits(end)*2.5e-2;
compLineArgs = {'Color','k','Linewidth',1};

for i = 1:2
    axx(i+2) = subplot(nRow,nCol, i + 2); hold on;
    ax = axx(i+2);
    switch i
        case 1
            val = TPR;
        case 2
            val = FPR;
    end
    % jitter points a little along the x axis
    plot([1+markerLineSpacing; 2-markerLineSpacing] + jitter, val', 'Color', col(3,:),'LineWidth',lwd_sub);
    plot(1 + jitter, val(:,1), 'o','MarkerSize',mks,'Color',col(1,:));
    plot(2 + jitter, val(:,2), 'o','MarkerSize',mks,'Color',col(2,:));
    
    ax.XAxis.Limits = [1,2] + 2*[-maxJitter,maxJitter];
    ax.XAxis.TickValues = [1, 2];
    ax.XAxis.TickLabels = instrumentOrder;
    
    ax.YAxis.Limits = [0,yLimits(i)];
    
    ax.YAxis.Label.String = labels{i};
    %
    if pval_TF(i) <= alpha0
        pltools.addComparisonLine(ax,1:2,y0(i),dy,compLineArgs{:});
        txt = pltools.equalText('p',pval_TF(i),1);
        text(1.5,y0(i) + dy,txt,...
            'Units','data','HorizontalAlignment','center',...
            'VerticalAlignment','bottom',...
            'FontSize',round(0.8*fts))
    end
end

pltools.formatAxisLabels(axx, fts, lwd);
axx(2).XAxis.FontSize = fts;
axx(3).XAxis.FontSize = fts;
axx(4).XAxis.FontSize = fts;

arrayfun(@(i) pltools.topLeftLabel(axx(i),topLabels{i},fts),1:4);
%
%