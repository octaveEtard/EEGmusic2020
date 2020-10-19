function topFolder = getPath(varargin)
%
% JoNmusic2020.linearBackwardModel
% Part of the JoNmusic2020 code.
% Author: Octave Etard
%
% This functions defines the base folder where the data to load, or results
% to save are located
%
switch varargin{1}
    
    case 'linearModelsResults'
        % where to save results, or load them from to plot
        topFolder = fullfile('..','..','data','linearModelsResults');
        
    case 'EEG'
        % where the EEG data is located
        % varargin{2}: 'raw' or 'processed'
        topFolder = fullfile('..','..','data','EEG',varargin{2});
        
    case 'features'
        % where the stimulus features are located
        topFolder = fullfile('..','..','data','stimuli','features');
        
    case 'MIDI'
        % where the MIDI files are located
        topFolder = fullfile('..','..','data','stimuli','rawInventions');
        
    case 'behav'
        % where the stimulus features are located
        topFolder = fullfile('..','..','data','behav');
        
    otherwise
        error('Unrecognised option.')
end
end
%
%