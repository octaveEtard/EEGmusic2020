function eegFolder = makePathEEGFolder(EEGopt)
%
% EEGmusic2020.makePathEEGFolder
% Part of the EEGmusic2020 code (github.com/octaveEtard/EEGmusic2020)
% Author: Octave Etard, 2020
%
% Generate path to EEG folder according to the following convention:
eegFolder = fullfile(EEGopt.baseFolder,sprintf('Fs-%i',EEGopt.Fs),EEGopt.proc,EEGopt.SID);
% e.g. someFolder/Fs-5000/HP-115/EBIP01
end
%
%