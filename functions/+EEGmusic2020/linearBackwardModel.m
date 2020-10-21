function CC = linearBackwardModel(conditions,SID,parts,EEGopt,...
    featureOpt,fields,opt,trainOpt)
%
% EEGmusic2020.linearBackwardModel
% Part of the EEGmusic2020 code.
% Author: Octave Etard
%
% Train backward models with a leave-one-data-part-out cross-validation and
% either leave-one-subject-out (generic models) or subject specific models.
%
% Do training + validation ; the validation performance can then be used to
% choose a model to test on independent held out data (not done inside this
% function), or a given model (corresponding to an a priori chosen
% regularisation coefficient) can be obtained and evaluated with this
% function.
%
nParts = numel(parts);
nSub = numel(SID);
Fs = EEGopt.Fs; % sampling rate

% --- options for the LMpackage functions
opt.getResponse = @EEGmusic2020.loadEEG; % function to load EEG ...
opt.getStimulus = @EEGmusic2020.loadFeature; % ... and stimulus feature
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
EEGopt.proc = EEGopt.procTrain;
featureFiles = EEGmusic2020.makePathFeatureFiles(conditions,parts,featureOpt,fields);
EEGFiles = EEGmusic2020.makePathEEGFiles(conditions,SID,parts,EEGopt);
% load data & compute compute cross matrices
[XtX,Xty,mX,mY,N] = LM.crossMatrices(featureFiles,EEGFiles,opt,'backward');

nCond = size(featureFiles,2);
% to store performance metrics
CC = cell(nParts,nSub,nCond);


%% --- Cross-validation
EEGopt.proc = EEGopt.procTest;

T_ = 0;
for iTestPart = 1:nParts
    if opt.printProgress
        fprintf('\nTest part = %i / %i\n',iTestPart,nParts);
        t1 = tic;
    end
    % iTestPart left out in the training data
    iTrainParts = [1:(iTestPart-1),(iTestPart+1):nParts];
    
    % using only training parts
    XtX_ = squeeze(sum(XtX(:,:,iTrainParts,:),3));
    % forming Xty_ on training data
    Xty_ = sum(Xty(:,:,iTrainParts,:),3);
    
    % Fitting model without mean : learning the mean on the training data
    % to remove it on the testing data
    mY_train = sum(N(iTrainParts) .* mY(:,iTrainParts),2) ./ sum(N(iTrainParts));
    
    for iTestSub = 1:nSub
        if opt.printProgress
            fprintf('\nTest sub = %i / %i\n',iTestSub,nSub);
        end
        if opt.generic % generic decoder
            % exclude iTestSub from training
            idxTrainSub = [1:(iTestSub-1),(iTestSub+1):nSub];
            
        else % subject-specific decoder
            idxTrainSub = iTestSub;
        end
        
        % sum over all train subjects
        XtX_train = sum(XtX_(:,:,idxTrainSub),3);
        Xty_train = sum(Xty_(:,:,:,idxTrainSub),4);
        
        % model fitted using all data bar iTestPart for all training
        % subjects
        tmp_model = LM.fitLinearModel(XtX_train,Xty_train,trainOpt);
        tmp_model = tmp_model.coeffs;
        
        mX_train = sum(N(iTrainParts) .* mX(:,iTrainParts,iTestSub),2) ./ sum(N(iTrainParts));
        
        % evaluating the model
        featureFiles_test = EEGmusic2020.makePathFeatureFiles(conditions,parts(iTestPart),featureOpt,fields);
        EEGFiles_test = EEGmusic2020.makePathEEGFiles(conditions,SID(iTestSub),parts(iTestPart),EEGopt);
        
        CC(iTestPart,iTestSub,:) = LM.testModel(tmp_model,featureFiles_test,EEGFiles_test,opt,'backward',mX_train,mY_train);
    end
    if opt.printProgress
        T_ = T_ + toc(t1);
        T = (nParts-iTestPart) * T_ / iTestPart;
        fprintf('\nEstimated remaining time = %s\n',seconds2HMS(T,{'h','mn','s'}));
    end
end
end
%
%