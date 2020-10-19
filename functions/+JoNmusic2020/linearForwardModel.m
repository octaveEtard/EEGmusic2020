function [model,CC] = linearForwardModel(conditions,SID,parts,EEGopt,...
    featureOpt,fields,opt,trainOpt)
%
% JoNmusic2020.linearForwardModel
% Part of the JoNmusic2020 code.
% Author: Octave Etard
%
% Train forward models with a leave-one-data-part-out cross-validation and
% either leave-one-subject-out (generic / population average models) or
% subject specific models.
%
% Do training + validation ; the validation performance can then be used to
% choose a model to test on independent held out data (not done inside this
% function), or a given model (corresponding to an a priori chosen
% regularisation coefficient) can be obtained and evaluated with this
% function.
%
%%
nSub = numel(SID);
nParts = numel(parts);
nCond = numel(conditions);
Fs = EEGopt.Fs; % sampling rate

% --- options for the LMpackage functions
opt.getResponse = @JoNmusic2020.loadEEG; % function to load EEG ...
opt.getStimulus = @JoNmusic2020.loadFeature; % ... and stimulus feature
% get stimulus will return a (nParts x nCond) variable; pool data from dim
% 2 (pool data across conditions (no effect if only 1 condition))
opt.sumFrom = 2;
% time to index
opt.minLag = floor(opt.minLagT * Fs);
opt.maxLag = ceil(opt.maxLagT * Fs);

opt.removeMean = true; % fitting a model without offset
opt.nPntsPerf = ceil(opt.perfSliceT * Fs) + 1; % time to index
opt.nStimPerFile = 1; % EEG dataset contains responses to 1 stim per file
opt.nChan = 2; % 2 channels in the EEG data
opt.sumStim = false; % does not matter as 1 stim per file
opt.unpad.do = false;
opt.predBatchSize = 'all';
opt.sumSub = false; % get data for all subjects to do cross-validation


%% --- form the XtX for all the data, as well as Xty for all subjects
% generate file paths
featureOpt.proc = featureOpt.procTrain;
featureFiles = JoNmusic2020.makePathFeatureFiles(conditions,parts,featureOpt,fields);
EEGFiles = JoNmusic2020.makePathEEGFiles(conditions,SID,parts,EEGopt);
% load data & compute compute cross matrices for all subjects
[XtX,Xty,mX,mY,N] = LM.crossMatrices(featureFiles,EEGFiles,opt,'forward');

% to store performance metrics
CC = cell(nParts,nSub,nCond);


%% --- Cross-validation
featureOpt.proc = featureOpt.procTest;

T_ = 0;
for iTestPart = 1:nParts
    if opt.printProgress
        fprintf('\nTest part = %i / %i\n',iTestPart,nParts);
        t1 = tic;
    end
    % iTestPart left out in the training data
    iTrainParts = [1:(iTestPart-1),(iTestPart+1):nParts];
    
    % using only training data
    XtX_train = (nSub-1) * sum(XtX(:,:,iTrainParts),3);
    % forming Xty_ on training data & averaging over channels
    Xty_ = squeeze(sum(Xty(:,:,iTrainParts,:),2:3) / opt.nChan);
    
    % Fitting model without mean : learning the mean on the training data
    % to remove it on the testing data
    mX_train = sum(N(iTrainParts) .* mX(:,iTrainParts),2) / sum(N(iTrainParts));
    
    for iTestSub = 1:nSub
        if opt.printProgress
            fprintf('\nTest sub = %i / %i\n',iTestSub,nSub);
        end
        if opt.generic % generic / population average models
            % exclude iTestSub from training
            idxTrainSub = [1:(iTestSub-1),(iTestSub+1):nSub];
            
        else % subject-specific models
            idxTrainSub = iTestSub;
        end
        
        % sum over all train subjects
        Xty_train = sum(Xty_(:,idxTrainSub),2);
        
        % model fitted using all data bar iTestPart for all training
        % subjects
        tmp_model = LM.fitLinearModel(XtX_train,Xty_train,trainOpt);
        tmp_model = tmp_model.coeffs;
        
        mY_train = sum(N(iTrainParts) .* mY(:,iTrainParts,idxTrainSub),1:3) ./ (opt.nChan * numel(idxTrainSub)*sum(N(iTrainParts)));
        
        % evaluating the model
        featureFiles_test = JoNmusic2020.makePathFeatureFiles(conditions,parts(iTestPart),featureOpt,fields);
        EEGFiles_test = JoNmusic2020.makePathEEGFiles(conditions,SID(iTestSub),parts(iTestPart),EEGopt);
        
        CC(iTestPart,iTestSub,:) = LM.testModel(tmp_model,featureFiles_test,EEGFiles_test,opt,'forward',mX_train,mY_train);
    end
    if opt.printProgress
        T_ = T_ + toc(t1);
        T = (nParts-iTestPart) * T_ / iTestPart;
        fprintf('\nEstimated remaining time = %s\n',seconds2HMS(T,{'h','mn','s'}));
    end
end


%%
% finally train will all data for each subject independently
% generic model = average of subject specific models at same lambda
% (since same stimulus)
XtX = sum(XtX,3);
Xty = squeeze(sum(Xty,2:3) / opt.nChan);

% fit model for all subjects at once: avoid having to inverse XtX multiple
% times (same for all subjects!): we treat the data for each subject as an
% independent output
model = LM.fitLinearModel(XtX,Xty,trainOpt);
model = model.coeffs;
end
%
%