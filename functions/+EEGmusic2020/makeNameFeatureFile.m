function featureFileName = makeNameFeatureFile( featureOpt )
%
% EEGmusic2020.makeNameFeatureFile
% Part of the EEGmusic2020 code.
% Author: Octave Etard, 2020
%
% Generate the name of the file containing the requested feature (file name
% with extension).
%
Fs = featureOpt.Fs; % sampling rate of the feature
proc = featureOpt.proc; % processing (filtering) of the feature
typeName = featureOpt.typeName; % feature name
condition = featureOpt.condition; % condition associated with this feature
iPart = featureOpt.part; % index of the part

% could be added to featureOpt if dynamic changes are necessary
vibType = 'zv'; % vibratos replaced by zeros
mixing = 1.25; % piano / guitar mixing in the competing conditions

switch condition
    case {'fGs','fPs'}
        instru = condition(2);
        audioName = sprintf('invent_%s_%s_%i',vibType,instru,iPart);
        
    % same mixed stimulus for fGc & fPc, but 2 feature files with
    % different name and attended / ignored field flipped
    case {'fGc','fPc'}
        audioName = sprintf('invent_%s_%s_PGf%s_%i',vibType,sprintf('%.2f',mixing),condition(2),iPart);  
end

% e.g. 'feature_Fs-5000-LP-2000-waveform_invent_zv_1.25_PGfP_2
featureFileName = sprintf('feature_Fs-%i-%s-%s_%s.mat',Fs,proc,typeName,audioName);

end
%
%