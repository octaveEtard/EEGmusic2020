function pathEEGFiles = makePathEEGFiles(conditions,allSID,parts,EEGopt)
%
% EEGmusic2020.makePathEEGFiles
% Part of the EEGmusic2020 code.
% Author: Octave Etard, 2020
%
% Generate paths to EEG files. This function returns a cell
% array of size (nParts x nConditions x nSubjects) containing paths to the
% files storing the relevant EEG data. 
%
% Each EEG file is associated with a (part x condition x subject)
% combination corresponding to the subject & stimulus that was played. Each
% file is moreover characterised by the processing that was applied to it
% (e.g. filering).
%
% EEGproc should either have the same number of elements as conditions, in
% this case it describes the EEG processing for each condition, or only one
% element which in this case applies to all conditions.
%
EEGproc = EEGopt.proc;
[conditions,allSID,EEGproc] = putInCell(conditions,allSID,EEGproc);  % FIXME

nCond = numel(conditions);
nParts = numel(parts);
nSub = numel(allSID);

if (1 < nCond) && (numel(EEGproc) == 1)
    EEGproc = repmat(EEGproc,1,nCond);
end

assert( numel(EEGproc) == nCond,...
    'The dimensions of ''EEGproc'' and ''conditions'' do not match.');


pathEEGFiles = cell(nParts,nCond,nSub);

for iCond = 1:nCond
    condition = conditions{iCond};
    
    for iiPart = 1:nParts
        iPart = parts(iiPart);
        
        for iSub = 1:nSub
            SID = allSID{iSub};
            % e.g. someFolder/Fs-5000/HP-115/EBIP01
            eegFolder = EEGmusic2020.makePathEEGFolder(EEGopt.baseFolder,EEGproc{iCond},EEGopt.Fs,SID);
            % e.g. 'HP-115-Fs-5000-EBIP01_fGs_2.set'
            eegFileName = EEGmusic2020.makeNameEEGDataFile(EEGproc{iCond},EEGopt.Fs,SID,condition,iPart,'.set');
            pathEEGFiles{iiPart,iCond,iSub} = {eegFolder,eegFileName};
        end
    end
end

pathEEGFiles = squeeze(pathEEGFiles);

if ismatrix(pathEEGFiles) && size(pathEEGFiles,1) == 1
    pathEEGFiles = pathEEGFiles';
end
end
%
%