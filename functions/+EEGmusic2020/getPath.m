function topFolder = getPath(varargin)
%
% EEGmusic2020.linearBackwardModel
% Part of the EEGmusic2020 code.
% Author: Octave Etard
%
% This functions defines the base folder where the data to load, or results
% to save are located.
%
% Currently assuming the following architecture:
%
% someFolder
% |--- EEGmusic2020.git                 (arbitray name)
% |    |--- someFolder
% |         |--- main analysis scripts
% |--- data                             (location of 'dataFolder')
%      |--- EEG
%      |    |--- processed
%      |--- stimuli
%      |    |--- features
%      |    |--- rawInventions
%      |--- behav
%
% location of the base folder containing the data ; % change to 
% 'path/to/top/folder/containing/data' if a different location is desired.
dataFolder = fullfile('..','..','data');
%
% See 'EEGmusic2020.makePathEEGFolder' 'EEGmusic2020.makePathEEGFolder'
% EEGmusic2020.'makePathFeatureFiles' EEGmusic2020.'makePathSaveResults'
% for the functions specifying the expected organisation of the data folder
%
% See 'EEGmusic2020.makeNameEEGDataFile' 'EEGmusic2020.makeNameEEGDataFile'
% for the functions specifying the expected file names.
%
switch varargin{1}
    
    case 'linearModelResults'
        % where to save results, or load them from to plot
        topFolder = fullfile(dataFolder,'linearModesResults');
        
    case 'EEG'
        % where the EEG data is located
        % varargin{2}: 'raw' or 'processed'
        topFolder = fullfile(dataFolder,'EEG',varargin{2});
        
    case 'features'
        % where the stimulus features are located
        topFolder = fullfile(dataFolder,'stimuli','features');
        
    case 'MIDI'
        % where the MIDI files are located
        topFolder = fullfile(dataFolder,'stimuli','rawInventions');
        
    case 'behav'
        % where the stimulus features are located
        topFolder = fullfile(dataFolder,'behav');
        
    otherwise
        error('Unrecognised option.')
end
end
%
%