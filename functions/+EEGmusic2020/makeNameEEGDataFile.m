function eegFileName = makeNameEEGDataFile(EEGopt)
%
% EEGmusic2020.makeNameEEGDataFile
% Part of the EEGmusic2020 code (github.com/octaveEtard/EEGmusic2020)
% Author: Octave Etard, 2020
%
% Generate name of EEG dataset according to the following convention:
eegFileName = sprintf('%s-Fs-%i-%s_%s_%i%s',EEGopt.proc,EEGopt.Fs,EEGopt.SID,EEGopt.condition,EEGopt.part,EEGopt.ext);
% e.g. 'HP-115-Fs-5000-EBIP01_fGs_2.set' 
end
%
%