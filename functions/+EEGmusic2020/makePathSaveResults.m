function [saveName,saveFolder] = makePathSaveResults(condition,EEGproc,...
    featureProc,featureTypeName,Fs,minLagT,maxLagT,modelType,baseSaveFolder)
%
% EEGmusic2020.makePathSaveResults
% Part of the EEGmusic2020 code.
% Author: Octave Etard, 2020
%
% Make file name and path to save results from forward or backward models.
%
[condition,EEGproc] = formatInput(condition,EEGproc);

saveName = sprintf('%s_%s_%s_%i_%ims.mat',modelType,condition,featureTypeName,1e3*minLagT,1e3*maxLagT);

saveFolder = fullfile(baseSaveFolder,...
    sprintf('Fs-%i',Fs),...
    EEGproc,...
    featureProc);
end

function varargout = formatInput(varargin)
for iarg = 1:nargin
    if iscell(varargin{iarg})
        f = sprintf('%s-',varargin{iarg}{:});
        varargin{iarg} = f(1:(end-1));
    end
end
varargout = varargin;
end
%
%