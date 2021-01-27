%
% figure_4
% Part of the EEGmusic2020 code (github.com/octaveEtard/EEGmusic2020)
% Author: Octave Etard, 2020
%
% Plot Figure 4.
%
% 1 - Load results of backward model from different time windows &
%       with mismatched stimulus / EEG
% 2 - Statistical analysis
% 2 - Plot
%
%% Load backward model results
conditions = {'fGs','fPs'};
featureTypeName = 'waveform';
EEGproc = 'HP-130';

featureProc = 'LP-2000';
Fs = 5000;

minLagT = [-15, 0, 15, 30, 0]*1e-3; % last one --> mismatched
maxLagT = minLagT+15*1e-3;

% analysing model with this regularisation coefficients
lambda0 = 10^-0.5;
nCond = numel(conditions);

nTime = numel(minLagT);

% top folder where the results are stored
baseSaveFolder = EEGmusic2020.getPath('linearModelResults');

% Load data from each time window & for each instrument
for it = 1:nTime
    for iCond = 1:nCond
        
        if it == nTime
            [saveName,saveFolder] = EEGmusic2020.makePathSaveResults(conditions{iCond},EEGproc,...
                featureProc,featureTypeName,Fs,minLagT(2),maxLagT(2),'mismatch_backward',baseSaveFolder);
        else
            [saveName,saveFolder] = EEGmusic2020.makePathSaveResults(conditions{iCond},EEGproc,...
                featureProc,featureTypeName,Fs,minLagT(it),maxLagT(it),'backward',baseSaveFolder);
        end
        
        d  = load(fullfile(saveFolder,saveName));
        
        if iCond == 1 && it == 1
            nSlices = sum(cellfun(@(m) size(m,1), d.CC(:,1)));
            nSub = numel(d.SID);
            nLambda = numel(d.train.method.lambda);
            
            CC = nan([nSlices,nSub,nLambda,nTime,nCond]);
            lambda = d.train.method.lambda;
        end
        % sanity check
        assert( all( d.train.method.lambda == lambda ) );
        tmp = vertcat(d.CC{:});
        CC(:,:,:,it,iCond) = reshape(tmp,[nSlices,nSub,nLambda]);
    end
end

iL0 = find(lambda == lambda0,1);
assert( ~isempty(iL0),'Cound not find requested lambda in data' );

% now slices x sub x time window x instruments
CC = squeeze( CC(:,:,iL0,:,:) );
subCC = squeeze(mean(CC,1)); % mean value for each subject

meanCC = squeeze(mean(subCC,1)); % mean value over the population
stdCC = squeeze(std(subCC,[],1)); % std over the population

nTime = nTime - 1; % last one corresponds to mismatched


%% --- Stats
% compare reconstruction for each subject in (0_15 ms) to (-15_0 ms)
pval_sub = nan(nSub,nCond);

for iCond = 1:nCond
    for iSub = 1:nSub
        % one-tailed paired Wilcoxon signed rank test
        pval_sub(iSub,iCond) = signrank(CC(:,iSub,2,iCond),CC(:,iSub,1,iCond),'tail','right');
    end
end

% compare reconstruction from (0_15 ms ; 15_30 ms ; 30_45ms) to (-15--0ms)
% at the population level
pval = nan(nTime,nCond);

for it = 2:nTime
    for iCond = 1:nCond
        % Wilcoxon paired signed rank test
        pval(it-1,iCond) = signrank(subCC(:,it,iCond),subCC(:,1,iCond),'tail','right');
    end
end

% compare reconstruction from (0_15 ms mismatched) to (-15--0ms) at the
% population level
pval_mism = nan(nCond,1);

for iCond = 1:nCond
    % Wilcoxon paired signed rank test
    pval_mism(iCond) = signrank(subCC(:,end,iCond),subCC(:,1,iCond),'tail','right');
end


% Subject level stats: FDR correction for multiple comparison over subjects
% x instruments
pval_sub_fdr = reshape(fdr(pval_sub(:)),size(pval_sub));

% Population level stats: FDR correction for multiple comparison over time
% windows x instruments
pval_fdr = reshape(fdr(pval(:)),size(pval));


%% Displaying some values in terminal
instru = {'Guitar','Piano'};
alpha0 = [0.1,1,5]/100; % significance levels
stars = {'***','**','*'}; % symbols corresponding to alpha0


fprintf('Subject level analysis:\n');
for iCond = 1:nCond
    n = sum(pval_sub_fdr(:,iCond) <= alpha0(end),1);
    fprintf('\t%s %i / %i subject(s) significant at %.3f\n',instru{iCond},n,nSub,alpha0(end));
end
n = pval_sub_fdr <= alpha0(end);
n = sum(n(:,1) & n(:,2),1);
fprintf('\tBoth instruments %i / %i subject(s) significant at %.3f\n\n',n,nSub,alpha0(end));

fprintf('Population level analysis:\n');
for iCond = 1:nCond
    for it = 2:nTime
        if alpha0(end) < pval_fdr(it-1,iCond)
            continue;
        end
        fprintf('\tWindow %i_%i ms for %s significant with p = %.2e\n',1e3*minLagT(it),1e3*maxLagT(it),instru{iCond},pval_fdr(it-1,iCond));
    end
end

fprintf('\nMismatched analysis:\n');
for iCond = 1:nCond
    if alpha0(end) < pval_mism(iCond)
        s = 'not ';
    else
        s = '';
    end
    fprintf('\t%s %i_%i ms mismatched %ssignificant with p = %.2e\n',instru{iCond},1e3*minLagT(end),1e3*maxLagT(end),s,pval_mism(iCond));
end

% order in which to plot subjects (bar plot)
[~,iSort] = sort(mean(subCC(:,2,:),3)); % mean correlation from 0_15 ms


%% Plotting
% --- for all plots
% fontsize and linewidth
fts = 11; lwd = 1;
% colours
[col,col_null] = EEGmusic2020.plotStyleArgs();
dyStars = 3e-3; % asterisks height above bars (data scale)

% --- bar plot plotting parameters
bw = 0.3; % histogram bar width (relative to space between subjects = 1)

% --- for both error bar plots
lwd_ = 1.5; % line width of error bars
% parameters of the line to indicate significance
compLineArgs = {'Color','k','Linewidth',lwd};

% --- error bar plot #1 param
errParam1 = struct();
% separation between null / meaningfull is 1 ; define other values relative
% to this
errParam1.instruStep = 0.7; % separation between instruments
errParam1.dx = 0.1; % subject line / errobar spacing
errParam1.xLimits = 0.25*[-1,1]; % white space beyond xMin/xMax
errParam1.compLineY = 6.5e-2; % comparison line height
% for lines linking subjects
errParam1.grey = col(3,:); % color
errParam1.lwd_sub = 0.6; % linewidth

% --- error bar plot #2 param
errParam2 = struct();
% separation between time windows is 1 ; define other values relative
% to this
errParam2.instruStep = 0.25; % separation between instruments
errParam2.xLimits = 0.15*[-1,1]; % white space beyond xMin/xMax
errParam2.compLineStep = 2e-2; % y spacing between comparison line
% first comparison line height (data)
errParam2.compLineY = errParam1.compLineY - errParam2.compLineStep;

% --- Y Axis limits: all plots
[m,M] = bounds(subCC,'all');
if 0 < m
    m = 0;
else
    m = 1.2*m;
end
M = max(1.2*M,M+3*dyStars);


% ---- Figure layout ---
axx = gobjects(5,1);

% width / height of each axis --> aspect ratio
xSpan = [7,7,2,3,2];
ySpan = [3,3,4,4,4];

% figure layout:
% <-  axx(1) ->
% <-  axx(2) ->
% axx(3) axx(4) axx(5)
nx = xSpan(1);
ny = sum(ySpan(1:3));

% --- actual plotting ---
fig = figure;
tl = tiledlayout(ny,nx);

% --- Bar plots
for iCond = 1:nCond
    
    axx(iCond) = nexttile(tl,[ySpan(iCond),xSpan(iCond)]); hold on;
    ax = axx(iCond);
    
    % 0_15 ms window
    bar((1:nSub)-bw/2,subCC(iSort,2,iCond),bw,...
        'EdgeColor','none','FaceColor',col(iCond,:));
    % -15_0 ms window (null)
    bar((1:nSub)+bw/2,subCC(iSort,1,iCond),bw,...
        'EdgeColor','none','FaceColor',col_null(iCond,:));
    
    % add asterisks depending on significance level
    for iSub = 1:nSub
        iSignifLevel = find(pval_sub_fdr(iSort(iSub),iCond) <= alpha0,1);
        
        if isempty(iSignifLevel)
            continue;
        end
        
        if 0 < subCC(iSort(iSub),2,iCond)
            y = subCC(iSort(iSub),2,iCond) + dyStars;
        else
            y = dyStars;
        end
        text(iSub,y,stars{iSignifLevel},'HorizontalAlignment','center',...
            'VerticalAlignment','middle',...
            'FontSize',fts);
    end
    
    ax.YAxis.Limits = [m,M];
    ax.XAxis.TickValues = (1:nSub);
    ax.XAxis.Limits = [1,nSub]+1.5*bw*[-1,1];
    
    text(1.5,M-4*dyStars,instru{iCond},'HorizontalAlignment','left',...
        'VerticalAlignment','middle',...
        'FontSize',fts);
end

ax.XAxis.Label.String = 'Subject';
ax.YAxis.Label.String = 'Correlation coefficient';
pltools.topLeftLabel(axx(1),'A',fts);

% --- Population comparison 0_15 ms vs -15_0 ms
axx(nCond+1) = nexttile(tl,[ySpan(nCond+1),xSpan(nCond+1)]); hold on;
ax = axx(nCond+1);

for iCond = 1:nCond
    
    x = (1:2) + (iCond-1) * (1+errParam1.instruStep);
    plot(x+errParam1.dx*[1,-1],subCC(:,1:2,iCond),'Color',errParam1.grey,...
        'LineWidth',errParam1.lwd_sub);
    
    % -15_0 ms window (null)
    errorbar(x(1),meanCC(1,iCond),stdCC(1,iCond),'o',...
        'Color',col_null(iCond,:),'LineWidth',lwd_,'MarkerFaceColor','w');
    % 0_15 ms window
    errorbar(x(2),meanCC(2,iCond),stdCC(2,iCond),'o',...
        'Color',col(iCond,:),'LineWidth',lwd_,'MarkerFaceColor','w');
    
    % add asterisks
    iSignifLevel = find(pval_fdr(1,iCond) <= alpha0,1);
    
    if isempty(iSignifLevel)
        continue;
    end
    p = stars{iSignifLevel};
    
    text(mean(x),errParam1.compLineY + dyStars,p,...
        'Units','data','HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'FontSize',fts);
    
    pltools.addComparisonLine(ax,x,errParam1.compLineY,dyStars,compLineArgs{:});
end

ax.XAxis.Limits = [1,x(end)]+errParam1.xLimits;
ax.XAxis.TickValues = [1.5,2.5+errParam1.instruStep];
ax.XAxis.TickLabels = instru;
ax.YAxis.Limits = [m,M];
ax.YAxis.Label.String = 'Correlation coefficient';
pltools.topLeftLabel(axx(nCond+1),'B',fts);

% --- Time window comparisons
axx(nCond+2) = nexttile(tl,[ySpan(nCond+2),xSpan(nCond+2)]); hold on;
ax = axx(nCond+2);

y0 = errParam2.compLineY;

for iCond = 1:nCond
    % -15_0 ms window (null)
    x = 1 + errParam2.instruStep*(iCond-1);
    errorbar(ax,x,meanCC(1,iCond),stdCC(1,iCond),'o',...
        'Color',col_null(iCond,:),'LineWidth',lwd_,'MarkerFaceColor','w');
    
    % Other windows
    x = (2:nTime) + errParam2.instruStep*(iCond-1);
    errorbar(ax,x,meanCC(2:nTime,iCond),stdCC(2:nTime,iCond),'o',...
        'Color',col(iCond,:),'LineWidth',lwd_,'MarkerFaceColor','w');
    
    % for easy legend
    lg(iCond) = plot(ax,x,meanCC(2:nTime,iCond),'o',...
        'Color',col(iCond,:),'LineWidth',lwd_,'MarkerFaceColor','w');
    
    % add pvalues
    x = (1:nTime) + errParam2.instruStep*(iCond-1);
    for iTime = 2:nTime
        if pval_fdr(iTime-1,iCond) <= alpha0(end)
            pltools.addComparisonLine(ax,x([1,iTime]),y0,dyStars,compLineArgs{:});
            
            p = pltools.equalText('p',pval_fdr(iTime-1,iCond),1);
            
            text(mean(x([1,iTime])),y0+dyStars,p,...
                'Units','data','HorizontalAlignment','center',...
                'VerticalAlignment','baseline',...
                'FontSize',round(0.8*fts));
            
            y0 = y0 + errParam2.compLineStep;
        end
    end
end

ax.XAxis.Limits = [1,x(end)]+errParam2.xLimits;
ax.XAxis.TickValues = (1:nTime) + errParam2.instruStep/2;
ax.XAxis.TickLabels = arrayfun(@(t1,t2) sprintf('%i to %i ms',1e3*t1,1e3*t2),...
    minLagT,maxLagT,'UniformOutput',false);
ax.XAxis.TickLabelRotation = 0;
ax.YAxis.Limits = [m,M];
legend(lg,instru,'box','off','Location','northeast','FontSize',fts);
pltools.topLeftLabel(axx(nCond+2),'C',fts);


% --- Population comparison 0_15 ms mismatched vs -15_0 ms
axx(nCond+3) = nexttile(tl,[ySpan(nCond+3),xSpan(nCond+3)]); hold on;
ax = axx(nCond+3);

for iCond = 1:nCond
    
    x = (1:2) + (iCond-1) * (1+errParam1.instruStep);
    plot(x+errParam1.dx*[1,-1],subCC(:,[1,end],iCond),'Color',errParam1.grey,...
        'LineWidth',errParam1.lwd_sub);
    
    % -15_0 ms window (null)
    errorbar(x(1),meanCC(1,iCond),stdCC(1,iCond),'o',...
        'Color',col_null(iCond,:),'LineWidth',lwd_,'MarkerFaceColor','w');
    % 0_15 ms mismatched window
    errorbar(x(2),meanCC(end,iCond),stdCC(end,iCond),'o',...
        'Color',col(iCond,:),'LineWidth',lwd_,'MarkerFaceColor','w');
    
    % add asterisks
    iSignifLevel = find(pval_mism(iCond) <= alpha0,1);
    
    if isempty(iSignifLevel)
        continue;
    end
    p = stars{iSignifLevel};
    
    text(mean(x),errParam1.compLineY + dyStars,p,...
        'Units','data','HorizontalAlignment','center',...
        'VerticalAlignment','middle',...
        'FontSize',fts);
    
    pltools.addComparisonLine(ax,x,errParam1.compLineY,dyStars,compLineArgs{:});
end

ax.XAxis.Limits = [1,x(end)]+errParam1.xLimits;
ax.XAxis.TickValues = [1.5,2.5+errParam1.instruStep];
ax.XAxis.TickLabels = instru;
ax.YAxis.Limits = [m,M];
pltools.topLeftLabel(axx(nCond+3),'D',fts);


pltools.formatAxisLabels(axx,fts,lwd);
tl.Padding = 'none';
tl.TileSpacing = 'none';

width = 15.25;
height = 10;
fileName = 'figure_4.pdf';
% pltools.printFigure(fig,'',fileName,600,width,height,1,1,1,'pdf',0);
%
%
