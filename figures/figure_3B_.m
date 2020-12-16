%
% figure_3B
% Part of the EEGmusic2020 code (github.com/octaveEtard/EEGmusic2020)
% Author: Octave Etard, 2020
%
% Plot Figure 3-B.
%
% 1 - Load forward model results
% 2 - Extract model at best CC
% 3 - Compare ROI vs baseline
% 4 - Plot
%
%% Load forward model results
condition = 'fGs-fPs';
featureTypeName = {'waveform','rectified'};
EEGproc = 'HP-80-HP-80';

featureProc = 'LP-2000';
Fs = 5000;

nChan = 2;
minLagT = -45e-3;
maxLagT = 100e-3;

nFeatures = numel(featureTypeName);

% top folder where the results are stored
baseSaveFolder = EEGmusic2020.getPath('linearModelResults');

for iFeature = 1:nFeatures
    
    [saveName,saveFolder] = EEGmusic2020.makePathSaveResults(condition,EEGproc,...
        featureProc,featureTypeName{iFeature},Fs,minLagT,maxLagT,'forward',baseSaveFolder);
    
    d  = load(fullfile(saveFolder,saveName));
    
    if iFeature == 1
        nSlices = sum(cellfun(@(m) size(m,1), d.CC(:,1)));
        nSub = numel(d.SID);
        lambda = d.train.method.lambda;
        nLambda = numel(lambda);
        
        CC = nan([nSlices,nSub,2,nChan,nLambda,nFeatures]);
        model = nan([size(d.model),nFeatures]);
    end
    
    tmp = vertcat(d.CC{:});
    CC(:,:,:,:,:,iFeature) = reshape(tmp,[nSlices,nSub,2,nChan,nLambda]);
    
    model(:,:,:,iFeature) = d.model;
end
clear d;


%% Analysis
% TRF time
minLag = floor(minLagT * Fs);
maxLag = ceil(maxLagT * Fs);
nLags = maxLag - minLag + 1;

tms = 1e3 * (-maxLag:1:-minLag) / Fs;

% plotting models at best regularisation coefficient
mCC = squeeze(mean(CC,1:4));
[maxCC,iMaxCC] = max(mCC,[],1);

figure; ax = axes(); hold on;
plot(lambda,mCC);
plot(lambda(iMaxCC),maxCC,'ko');
ax.XAxis.Scale = 'log';

% lambda0 = 10^1;
% [~,iL0] = min(abs(lambda-lambda0));
% % iMaxCC = iMaxCC(2) * [1,1];
% iMaxCC =  [iL0,iMaxCC(2)];
% plot(lambda(iMaxCC),maxCC,'go');


trf = nan(nLags,nSub,2);
for iFeature = 1:2
    trf(:,:,iFeature) = model(:,:,iMaxCC(iFeature),iFeature);
end

meanTRF = squeeze(mean(trf,2));

% baseline for significance
tmin = -maxLagT * 1e3;
tmax = 0;
bs = tmin <= tms & tms <= tmax;

nbs = sum(bs); % number of points in the baseline
ROI = ~bs; % region of interest
nROI = sum(ROI); % number of points in the region of interest

% testing ROI vs baseline
pval = nan(nROI,nFeatures);

mu = nan(2,1);
sigma = nan(2,1);

for iFeature = 1:2
    % get null distribution by fitting Gaussian pdf to the pooled values
    % from the baseline
    nd = meanTRF(bs,iFeature);
    nd = fitdist(nd,'Normal'); % null distribution
    
    % two-tailed uncorrected pval based on nd
    roi = meanTRF(ROI,iFeature);
    pval(:,iFeature) = 2*min(normcdf(roi,nd.mu,nd.sigma),1-normcdf(roi,nd.mu,nd.sigma));
end

% significance level
alpha0 = 5/100;
% FDR correction: corrected significance level (time points x instruments)
alpha = fdr( pval(:),alpha0 );

% plot shading at +/- 1 std around the mean
stdTRF = squeeze(std(trf,[],2));


%% Plotting
fts = 11; lwd = 1;

lgValues = {'TFS','Note onsets removed'};
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

ax.XAxis.Limits = [-100,45];
ax.YAxis.TickValues = -1:1;

pltools.formatAxisLabels(ax,fts,lwd);
pltools.topLeftLabel(ax,'B',fts);
%
%
