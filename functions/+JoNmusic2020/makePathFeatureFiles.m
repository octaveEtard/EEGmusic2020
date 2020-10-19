function pathFeatureFiles = makePathFeatureFiles(conditions,parts,featureOpt,fields)
%
% JoNmusic2020.makePathFeatureFiles
% Part of the JoNmusic2020 code.
% Author: Octave Etard, 2020
%
% Generate paths to stimulus feature files. This function returns a cell
% array of size (nParts x nConditions) containing paths to the files
% storing the relevant features. 
%
% A feature corresponds to a particular representaion of a stimulus
% (e.g. waveform, envelope). Each stimulus is associated with a
% (condition x part) combination corresponding to the EEG recording where
% it was played. Each stimulus can be described by several features. Each
% feature is described by a name (e.g. waveform, envelope) and a processing
% applied to it (e.g. high-pass filering).
%
% Here it is assumed that all stimuli will be described by the same number
% of features ('nFeatures'), and the same 'kind' of feature (name x
% processing). Each feature is associated with a file from which a number
% of 'fields' can be selected (e.g. 'attended', 'ignored'). Here the same
% fields will be selected for all features.
%
% In the end nFeatures x nFields variables will be prepared.
%

% Make sure input are cell arrays if not already the case. 
[conditions,featureName,featureProc,fields] = putInCell( conditions,...
    featureOpt.typeName,...
    featureOpt.proc,...
    fields);
% featureName / featureProc apply to all stimuli (condition x part
% combinations). They should be cell arrays with the same number of
% elements.
%
if ~iscell(fields)
    fields = {fields};
end
% fields apply to all features. It should be a cell array containing the
% fields to select from each feature file.
%
nFeatures = numel(featureName);
assert( nFeatures == numel(featureProc) );

% for each condition x part, nFeatures paths will be generated
nCond = numel(conditions);
nParts = numel(parts);

Fs = featureOpt.Fs;

pathFeatureFiles = cell(nParts,nCond);

for iCond = 1:nCond
    featureOpt.condition = conditions{iCond};
    
    for iiPart = 1:nParts
        iPart = parts(iiPart);
        
        c = cell(nFeatures,1);
        for iFeature = 1:nFeatures
            featureOpt.typeName = featureName{iFeature};
            featureOpt.proc = featureProc{iFeature};
            
            Fs_ = sprintf('Fs-%i',Fs);
            % e.g. someFolder/Fs-5000/LP-2000/waveform
            envFolder = fullfile(featureOpt.baseFolder,Fs_,featureOpt.proc,featureOpt.typeName);
            % e.g. 'feature_Fs-5000-LP-2000-waveform_invent_zv_1.25_PG_2
            envFileName = JoNmusic2020.makeNameFeature(featureOpt,iPart);
            
            c{iFeature} = {fullfile(envFolder,envFileName),fields};
        end
        pathFeatureFiles{iiPart,iCond} = c;
    end
end

pathFeatureFiles = squeeze(pathFeatureFiles);

if ismatrix(pathFeatureFiles) && size(pathFeatureFiles,1) == 1
    pathFeatureFiles = pathFeatureFiles';
end
end
%
%