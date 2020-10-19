function [EEG,iB,iE] = loadEEG(c)
%
% JoNmusic2020.loadEEG
% Part of the JoNmusic2020 code.
% Author: Octave Etard, 2020
%
% Load EEG data set, remove 'Sound' channel, and return EEG data matrix,
% iB index of point where stimulus starts, and iE where it ends.
%
% Requires the EEGLAB toolbox ('pop_loadset' function).
%
% Input:
%    c{1}{1}: folder containing file
%    c{1}{2}: file name
%
% Output:
%   EEG: data matrix nPnts x nChannels
%   iB: index in EEG where the stimulus starts
%   iE: index in EEG where the stimulus ends
%
c = c{1};
EEG = pop_loadset(c{2}, c{1});

iSound = findIndexSoundChannel(EEG);
% actual EEG channels;
iEEG = 1:EEG.nbchan;
iEEG(iSound) = [];

% index of point in EEG where the stimulus started
iB = findLatencyEvent(EEG,'stimBegin');
if isempty(iB)
    error('Could not find stimulus onset');
end
% index of point in EEG where the stimulus ended
iE = findLatencyEvent(EEG,'stimEnd');
if isempty(iB)
    error('Could not find stimulus end');
end

EEG = EEG.data(iEEG,:)';
end
%
%