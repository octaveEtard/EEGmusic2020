%
% music_SI_forward
% Part of the EEGmusic2020 code.
% Author: Octave Etard, 2020
%
% Train linear forward models in the Single Instrument conditions with
% leave-one-data-part-out and leave-one-subject-out cross-validation; then
% save the results.
%
%
%% Parameters
% subject IDs to use
SID = 1:17;
% 'EBIP01' to 'EBIP17'
SID = arrayfun(@(idx) sprintf('EBIP%02i',idx),SID,'UniformOutput',false);
% data from each condition divided in parts corresponding to different
% inventions: data parts indices to use here
% (1 was training block, not used)
parts = 2:7;
% conditions:
% shorthands for the SI conditions 'focus(Guitar/Piano)Single'
conditions = {'fGs','fPs'};
% sampling rate
Fs = 5000;
% processing of the EEG to use
EEGopt = [];
EEGopt.proc = 'HP-130'; % high-passed at 130 Hz

% name of the feature describing the stimulus
featureOpt = [];
featureOpt.typeName = 'waveform';
% processing of the feature
featureProc = 'LP-2000'; % low-passed at 2000 Hz (anti-aliasing / resampling)
fields = 'attended'; % only 1 instrument

% time window in which to train the model ; understood as time lag of
% predictor (here stimulus) with respect to predicted (here EEG) -->
% negative latencies = stimulus preceding EEG = causal / meaningful
opt = [];
opt.minLagT = -45e-3; % in s
opt.maxLagT = 100e-3;

% estimate performance using 10-s slices
opt.perfSliceT = 10; % in s

% top folder where to store the results
baseSaveFolder = EEGmusic2020.getPath('linearModelResults');


%%
nCond = numel(conditions);

% where is the EEG / stimulus feature data located
EEGopt.baseFolder = EEGmusic2020.getPath('EEG','processed');
featureOpt.baseFolder = EEGmusic2020.getPath('features');

EEGopt.Fs = Fs;
featureOpt.Fs = Fs;
featureOpt.procTrain = featureProc;
featureOpt.procTest = featureProc;

% --- options for the LMpackage functions
opt.nFeatures = 1; % 1 predictor feature

% ---
% train a generic / population average model
opt.generic = true;
% display some progress information
opt.printProgress = true;

% --- options for the ridge regression
trainOpt = [];
trainOpt.printOut = false;
trainOpt.accumulate = true;
trainOpt.method.name = 'ridge-eig-XtX'; % use ridge regression
trainOpt.method.lambda = 10.^(-6:0.5:6); % regularisation parameters
trainOpt.method.normaliseLambda = true;


%% For each condition train & test model
for iCond = 1:nCond
    if opt.printProgress
        fprintf('Starting cond %i / %i\n\n',iCond,nCond);
    end
    % train & test (cross-validation) a forward model for each intrument
    % / condition (guitar or piano)
    [model,CC] = EEGmusic2020.linearForwardModel(conditions(iCond),SID,...
        parts,EEGopt,featureOpt,fields,opt,trainOpt);
    
    %% save results
    d = [];
    % parameters
    d.SID = SID;
    d.condition = conditions{iCond};
    d.parts = parts;
    
    d.EEG = EEGopt;
    d.feature = featureOpt;
    d.fields = fields;
    d.opt = opt;
    d.train = trainOpt;
    
    % results
    d.CC = CC;
    d.model = model;
    
    [saveName,saveFolder] = EEGmusic2020.makePathSaveResults(conditions{iCond},EEGopt.proc,...
        featureProc,featureOpt.typeName,Fs,opt.minLagT,opt.maxLagT,'forward',baseSaveFolder);
    
    LM.save(d,saveName,saveFolder);
end
%
%