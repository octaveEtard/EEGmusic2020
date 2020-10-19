%
% music_SI_backward
% Part of the JoNmusic2020 code.
% Author: Octave Etard, 2020
%
% Train linear backward models in the Single Instrument conditions for each
% subject with leave-one-data-part-out cross-validation; then save the
% results.
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
% short hands for the SI conditions 'focus(Guitar/Piano)Single'
conditions = {'fGs','fPs'};
% sampling rate
Fs = 5000;
% processing of the EEG to use
EEGproc = {'HP-175','HP-115'}; % high-passed at: guitar: 175 / piano: 115 Hz

% name of the feature describing the stimulus
featureOpt = struct();
featureOpt.typeName = 'waveform';
% processing of the feature
featureOpt.proc = 'LP-2000'; % low-passed at 2000 Hz (anti-aliasing / resampling)
fields = 'attended'; % only 1 instrument

% time window in which to train the model ; understood as time lag of
% predictor (here EEG) with respect to predicted (here stimulus) -->
% positive latencies = EEG post stimulus = causal / meaningful
opt = struct();
opt.minLagT = 0;
opt.maxLagT = opt.minLagT + 15e-3; % 15 ms reconstruction window

% estimate performance using 10-s slices
opt.perfSliceT = 10; % in s

% top folder where to store the results
baseSaveFolder = JoNmusic2020.getPath('linearModelsResults');


%%
nCond = numel(conditions);

% where is the EEG / stimulus data located
EEGopt = struct();
EEGopt.baseFolder = JoNmusic2020.getPath('EEG','processed');
featureOpt.baseFolder = JoNmusic2020.getPath('features');

EEGopt.Fs = Fs;
featureOpt.Fs = Fs;

% --- options for the LMpackage functions
opt.nFeatures = 1; % 1 predicting 1 stimulus feature

% ---
% train subject specific decoders
opt.generic = false;
% display some progress information
opt.printProgress = true;

% --- options for the ridge regression 
trainOpt = struct();
trainOpt.printOut = false;
trainOpt.accumulate = true;
trainOpt.method.name = 'ridge-eig-XtX'; % ridge regression
trainOpt.method.lambda = 10^-0.5; % regularisation parameters
trainOpt.method.normaliseLambda = true;


%% For each condition train & test model
for iCond = 1:nCond
     if opt.printProgress
        fprintf('Starting cond %i / %i\n\n',iCond,nCond);
    end   
    % train & test (cross-validation) a backward model for each intrument
    % / condition  & subjects (guitar or piano)
    EEGopt.procTrain = EEGproc{iCond};
    EEGopt.procTest = EEGproc{iCond};

    CC = JoNmusic2020.linearBackwardModel(conditions{iCond},SID,parts,...
        EEGopt,featureOpt,fields,opt,trainOpt);

    %% save results
    d = struct();
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
    
    [saveName,saveFolder] = JoNmusic2020.makePathSaveResults(conditions{iCond},EEGproc{iCond},...
        featureOpt.proc,featureOpt.typeName,Fs,opt.minLagT,opt.maxLagT,'backward',baseSaveFolder);
    
    LM.save(d,saveName,saveFolder);
end
%
%