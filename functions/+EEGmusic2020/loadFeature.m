function features = loadFeature(c)
%
% EEGmusic2020.loadFeature
% Part of the EEGmusic2020 code.
% Author: Octave Etard, 2020
%
% Load feature files, extract relevant fields and create 'features' matrix.
%
% Input:
%   c: cell array such that:
%   c{1}: cell array of size nFeatures x 1, such for each element:
%   c{1}{iFeature}: cell array of size 1 x 2 with
%   c{1}{iFeature}{1}: path to feature file
%   c{1}{iFeature}{2}: fields to load from that file
%
% Ouput:
%   features: matrix of size nPnts x nCol with:
%       nPnts = number of points in all features
%       nCol = sum( iFeature * nFields(iFeature))
% [feature_1 <field_1..field_n1> feature_2 <field_1..field_n2> ...]
c = c{1};
nFeatures = numel(c);

% size of output
nCol = sum(cellfun(@(c) numel(c{2}),c));
iCol = 1;

for iFeature = 1:nFeatures
    
    filePath = c{iFeature}{1};
    fields = c{iFeature}{2};
    d = load(filePath);
    
    if iFeature == 1
        nPnts = size(d.(fields{1}),1); % should be the same for all features!
        features = zeros(nPnts,nCol);
    end
    nFields = numel(fields);
    for iUse = 1:nFields
        features(:,iCol) = d.(fields{iUse});
        iCol = iCol + 1;
    end
end
end
%
%