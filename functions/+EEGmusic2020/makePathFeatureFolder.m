function featureFolder = makePathFeatureFolder(featureOpt)
%
% EEGmusic2020.makePathFeatureFolder
% Part of the EEGmusic2020 code (github.com/octaveEtard/EEGmusic2020)
% Author: Octave Etard, 2020
%
% Generate path to feature folder according to the following convention:
featureFolder = fullfile(featureOpt.baseFolder,sprintf('Fs-%i',featureOpt.Fs),featureOpt.proc,featureOpt.typeName);
end
%