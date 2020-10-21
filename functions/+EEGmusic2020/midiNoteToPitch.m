function notes = midiNoteToPitch(notes,baseF)
%
% EEGmusic2020.midiNoteToPitch
% Part of the EEGmusic2020 code (github.com/octaveEtard/EEGmusic2020)
% Author: Octave Etard, 2020
%
% Get pitch in Hz of midi notes (convert MIDI note indices to pitch value)
if nargin < 2
    baseF = 440;
end
notes = baseF * 2.^((notes-69)/12);
end
%
%