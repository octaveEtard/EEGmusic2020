%
% figure_6
% Part of the EEGmusic2020 code (github.com/octaveEtard/EEGmusic2020)
% Author: Octave Etard, 2020
%
% Plot Figure 6.
%
% 1 - Load results from backward models.
% 2 - Load results from forward models.
% 3 - Plot.
%
%% Load backward model results
featureTypeName = 'waveform';
EEGproc = 'HP-130';
instrumentOrder = {'Guitar','Piano'};
attention = {'attended','ignored'};

Fs = 5000;
featureProc = 'LP-2000';

% analysing model with this regularisation coefficients
lambda0 = 10^-0.5;

minLagT = 0;
maxLagT = 15e-3;
nInstru = 2;

% top folder where the results are stored
baseSaveFolder = EEGmusic2020.getPath('linearModelResults');

for iInstru = 1:nInstru
    
    condition = sprintf('f%sc-AI',instrumentOrder{iInstru}(1));
    [saveName,saveFolder] = EEGmusic2020.makePathSaveResults(condition,EEGproc,...
        featureProc,featureTypeName,Fs,minLagT,maxLagT,'backward',baseSaveFolder);
    
    d  = load(fullfile(saveFolder,saveName));
    
    if iInstru == 1
        nSlices = sum(cellfun(@(m) size(m,1), d.CC(:,1)));
        nSub = numel(d.SID);
        lambda = d.train.method.lambda;
        nLambda = numel(lambda);
        % slices x nSub x nLambda x attention x condition
        CC = nan([nSlices,nSub,nLambda,2,nInstru]);
    end
    
    tmp = vertcat(d.CC{:});
    CC(:,:,:,:,iInstru) = reshape(tmp,[nSlices,nSub,nLambda,2]);
end
% now CC contains (attended / ignored) x condition, that is:
%    piano A / guitar I for fGc
%   guitar A / piano I for fPc
% we want (attended / ignored) x instrument instead, that is:
%    A / I for piano
%    A / I for guitar
% hence:
CC(:,:,:,2,:) = CC(:,:,:,2,[2,1]);

iL0 = find(lambda == lambda0,1);
assert(lambda(iL0) == lambda0,'Could not find requested lambda in data');

% mean CC for each subject at desired lambda
subCC = squeeze(mean(CC(:,:,iL0,:,:),1));
% population mean & std
meanCC = squeeze(mean(subCC,1));
stdCC = squeeze(std(subCC,[],1));


%% --- Stats
% attention test
pval_backward = nan(nInstru,1);
for iInstru = 1:nInstru
    % Wilcoxon paired signed rank test for zero median.
    pval_backward(iInstru) = signrank(subCC(:,1,iInstru),subCC(:,2,iInstru));
end
% % FDR correction
pval_backward_fdr = fdr(pval_backward)
% pval_backward_fdr = fdr_th(pval_backward)


%% Load forward model results
condition = 'fGc-fPc'; % pooled data over the 2 CI conditions
EEGproc = 'HP-130';
featureProc = 'LP-2000';
minLagT = -45e-3;
maxLagT = 100e-3;


[saveName,saveFolder] = EEGmusic2020.makePathSaveResults(condition,EEGproc,...
    featureProc,featureTypeName,Fs,minLagT,maxLagT,'forward',baseSaveFolder);

d = load(fullfile(saveFolder,saveName));

nSlices = sum(cellfun(@(m) size(m,1), d.CC(:,1)));
nChan = 2;
nSub = numel(d.SID);
lambda = d.train.method.lambda;
nLambda = numel(lambda);

CC_forward = reshape(vertcat(d.CC{:}),[nSlices,nSub,2,nChan,nLambda]);

minLag = floor(minLagT * Fs);
maxLag = ceil(maxLagT * Fs);
nLags = maxLag - minLag + 1;
tms = 1e3 * (-maxLag:1:-minLag) / Fs;

% plotting model at best regularisation coefficient
mCC_forward = squeeze(mean(CC_forward,1:4));
[~,iMaxCC_forward] = max(mCC_forward,[],1);

trf = reshape(d.model,[nLags,2,nSub,nLambda]);
% time x subjects x attention
trf = permute(trf(:,:,:,iMaxCC_forward),[1,3,2]);

% where to compare attended & ignored TRFs
ROI = 0 < tms & tms <= 15; % region of interest
nROI = sum(ROI); % number of points in the region of interest
idxROI = find(ROI); % index of points in ROI

pval_forward = nan(nROI,1);
% attention test
for it = 1:nROI
    pval_forward(it) = signrank(trf(idxROI(it),:,1),trf(idxROI(it),:,2));
end
% significance level
alpha0 = 1/100;
% FDR correction: corrected significance level
alpha_forward = fdr( pval_forward, alpha0 );

% average over subjects
meanTRF = squeeze(mean(trf,2));
% plot shading at +/- 1 std around the mean
stdTRF = squeeze(std(trf,[],2));


%% Plotting
% --- for all plots
% fontsize and linewidth
fts = 11; lwd = 1;
% default colour & plotting parameters
[col,col_null,shadingArgs,nsArgs,sgArgs] = EEGmusic2020.plotStyleArgs();

% --- bawckward model plotting parameters
LBMParam = struct();

% subject lines
LBMParam.sub_lwd = 0.6; % width
LBMParam.sub_color = col(3,:); % color
LBMParam.sub_dx = 0.1; % horizontal offset
LBMParam.eb_lwd = 1.5; % error bar linewidth
LBMParam.yText = 4.5e-2; % text height


axx = gobjects(2,1);

fig = figure;
tl = tiledlayout(1,2);

% --- Panel A (backward models)
axx(1) = nexttile(tl); hold on;
ax = axx(1);

for iInstru = 1:nInstru
    
    x = (1:2) + (iInstru-1) * 2;
    plot(x+LBMParam.sub_dx*[1,-1],subCC(:,:,iInstru),...
        'Color',LBMParam.sub_color,...
        'LineWidth',LBMParam.sub_lwd );
    errorbar(x,meanCC(:,iInstru),stdCC(:,iInstru),'o',...
        'Color',col(iInstru,:),...
        'LineWidth',LBMParam.eb_lwd ,...
        'MarkerFaceColor','w');
    
    text(mean(x),LBMParam.yText,instrumentOrder{iInstru},...
        'Units','data','HorizontalAlignment','center',...
        'VerticalAlignment','bottom',...
        'FontSize',fts);
end

ax.XAxis.Limits = [0.5,4.5];
ax.YAxis.Limits = [-2,5]/100;
ax.YAxis.Label.String = 'Correlation coefficient';

ax.XAxis.TickValues = 1:4;
ax.XAxis.TickLabels = {'Attended','Ignored','Attended','Ignored'};


% --- Panel B (forward models)
lsty = {'-','-'};
lgValues = {'Attended','Ignored'}; % legend

signif = pval_forward <= alpha_forward;

% for display purposes
m = max(abs([meanTRF+stdTRF,meanTRF-stdTRF]),[],'all');
meanTRF = meanTRF / m;
stdTRF = stdTRF / m;

axx(2) = nexttile(tl); hold on;
ax = axx(2);

[sgHandles,nsHandles] = EEGmusic2020.plot_trf(ax,tms,meanTRF,stdTRF,ROI,signif,...
    {},nsArgs,shadingArgs,...
    col,col,col_null,lsty,fts,lgValues);

ax.XAxis.Limits = [-10,45];
ax.YAxis.TickValues = [-1,1];

pltools.formatAxisLabels(axx,fts,lwd);
tl.Padding = 'none';
tl.TileSpacing = 'none';

pltools.topLeftLabel(axx(1),'A',fts);
pltools.topLeftLabel(axx(2),'B',fts);

width = 15.25;
height = 5;
fileName = 'figure_6.pdf';
% pltools.printFigure(fig,'',fileName,600,width,height,1,1,1,'pdf',0);
%
%