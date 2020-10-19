function col = getColorOrder(addBlack,nCol)
%
% pltools.printFigure
% Part of the pltools package (github.com/octaveEtard/pltools).
% Author: Octave Etard, 2020
%
if nargin < 1
    addBlack = 0;
end
if nargin < 2
    nCol = 7;
end
% 0: do no add black
% i: add black at position 'i'

if nargin < 2
    if addBlack > 0
        nCol = 8;
    else
        nCol = 7;
    end
end

col = [0, 0.4470, 0.7410;...
    0.8500, 0.3250, 0.0980;...
    0.9290, 0.6940, 0.1250;...
    0.4940, 0.1840, 0.5560;...
    0.4660, 0.6740, 0.1880;...
    0.3010, 0.7450, 0.9330;...
    0.6350, 0.0780, 0.1840];

addBlack = min(addBlack,size(col,1)+1);

if addBlack > 0
    col = [col(1:(addBlack-1),:);...
        0,0,0;...
        col(addBlack:end,:)];
end

nRep = ceil( nCol / size(col,1) );

if nRep > 1
    col = repmat(col,nRep,1);
end
    
col = col(1:nCol,:);

end