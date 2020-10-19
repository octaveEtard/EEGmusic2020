function eegFileName = makeNameEEGDataFile(proc,Fs,SID,condition,iPart,ext)
%
% JoNmusic2020.makeNameEEGDataFile
% Part of the JoNmusic2020 code (github.com/octaveEtard/JoNmusic2020)
% Author: Octave Etard, 2020
%
% Generate name of EEG dataset according to the following convention:
eegFileName = sprintf('%s-Fs-%i-%s_%s_%i%s',proc,Fs,SID,condition,iPart,ext);
% e.g. 'HP-115-Fs-5000-EBIP01_fGs_2.set' 
end
%
%