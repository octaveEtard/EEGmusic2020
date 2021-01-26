%
% figure_3A
% Part of the EEGmusic2020 code (github.com/octaveEtard/EEGmusic2020)
% Author: Octave Etard, 2020
%
% Plot Figure 3-A.
%
% 1 - Load forward model results
% 2 - Extract model at best CC
% 3 - Compare ROI vs baseline
% 4 - Plot
%
%% Load forward model results
conditions = {'fGs','fPs'};
featureTypeName = 'waveform';
EEGproc = 'HP-130';

featureProc = 'LP-2000';
Fs = 5000;

nChan = 2;
minLagT = -45e-3;
maxLagT = 100e-3;

nCond = numel(conditions);

% top folder where the results are stored
baseSaveFolder = EEGmusic2020.getPath('linearModelResults');

for iCond = 1:nCond
    
    [saveName,saveFolder] = EEGmusic2020.makePathSaveResults(conditions{iCond},EEGproc,...
    featureProc,featureTypeName,Fs,minLagT,maxLagT,'forward',baseSaveFolder);
    
    d  = load(fullfile(saveFolder,saveName));
    
    if iCond == 1
        nSlices = sum(cellfun(@(m) size(m,1), d.CC(:,1)));
        nSub = numel(d.SID);
        lambda = d.train.method.lambda;
        nLambda = numel(lambda);
        
        CC = nan([nSlices,nSub,nChan,nLambda,nCond]);
        model = nan([size(d.model),nCond]);
    end
    
    tmp = vertcat(d.CC{:});
    CC(:,:,:,:,iCond) = reshape(tmp,[nSlices,nSub,nChan,nLambda]);

    model(:,:,:,iCond) = d.model;
end

clear d;


%% Analysis
% TRF time
minLag = floor(minLagT * Fs);
maxLag = ceil(maxLagT * Fs);
nLags = maxLag - minLag + 1;

tms = 1e3 * (-maxLag:1:-minLag) / Fs;

% plotting models at best regularisation coefficient
mCC = squeeze(mean(CC,1:3));
[maxCC,iMaxCC] = max(mCC,[],1);

trf = nan(nLags,nSub,2);
for iInstru = 1:2
    trf(:,:,iInstru) = model(:,:,iMaxCC(iInstru),iInstru);
end
% average over subjects
meanTRF = squeeze(mean(trf,2));

% baseline for significance
tmin = -maxLagT * 1e3;
tmax = 0;
bs = tmin <= tms & tms <= tmax;

nbs = sum(bs); % number of points in the baseline
ROI = ~bs; % region of interest
nROI = sum(ROI); % number of points in the region of interest

% testing ROI vs baseline
pval = nan(nROI,nCond);

mu = nan(2,1);
sigma = nan(2,1);

for iInstru = 1:2
    % get null distribution by fitting Gaussian pdf to the pooled values
    % from the baseline
    nd = meanTRF(bs,iInstru);
    nd = fitdist(nd,'Normal'); % null distribution
    
    % two-tailed uncorrected pval based on nd
    roi = meanTRF(ROI,iInstru);
    pval(:,iInstru) = 2*min(normcdf(roi,nd.mu,nd.sigma),1-normcdf(roi,nd.mu,nd.sigma));
end

% significance level
alpha0 = 5/100;
% FDR correction: corrected significance level (time points x instruments)
alpha = fdr( pval(:),alpha0 );

% plot shading at +/- 1 std around the mean
stdTRF = squeeze(std(trf,[],2));


%% Plotting
fts = 11; lwd = 1;

lgValues = {'Guitar','Piano'};
lgHandle = gobjects(2,1);

% colors and linestyle
[col,col_null,shadingArgs,nsArgs,sgArgs] = EEGmusic2020.plotStyleArgs();
lsty = {'-','-.'};

signif = pval <= alpha;

yScaleFactor = max(abs(meanTRF)); % for display purposes
meanTRF = meanTRF ./ yScaleFactor;
stdTRF = stdTRF ./ yScaleFactor;

fig = figure;
ax = axes(); hold on;

[sgHandles,nsHandles] = EEGmusic2020.plot_trf(ax,tms,meanTRF,stdTRF,ROI,signif,...
    sgArgs,nsArgs,shadingArgs,...
    col,col_null,col_null,lsty,fts,lgValues);

% making sure dashed-red is on top of continuous-black
uistack(nsHandles([2,1]),'top');
uistack(sgHandles([2,1]),'top');

ax.XAxis.Limits = [-10,45];
ax.XAxis.TickValues = -10:10:40;
ax.YAxis.TickValues = -1:1;

pltools.formatAxisLabels(ax,fts,lwd);
pltools.topLeftLabel(ax,'A',fts);

width = 15.25/2;
height = 4;
fileName = 'figure_3A';
pltools.printFigure(fig,'',fileName,600,width,height,1,1,1,'pdf',1);
%
%