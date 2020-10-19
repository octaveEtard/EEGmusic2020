function p = shadeXarea(ax,y,xMin,xMax,varargin)
%
% pltools.shadeXarea
% Part of the pltools package (github.com/octaveEtard/pltools).
% Author: Octave Etard, 2020
%
% shade area alongside y, between xMin and xMax
assert(length(xMin) == length(xMax));

if length(y) ~= length(xMin)
    assert(length(xMin) == 1);
    y = [y(1),y(end)];
    xMin = [xMin,xMin];
    xMax = [xMax,xMax];
end
% row vector
y = y(:)';
xMin = xMin(:)';
xMax = xMax(:)';

y = [y,fliplr(y)];
x = [xMin(:)',fliplr(xMax(:)')];
 
p = patch(ax,'XData',x,'YData',y,varargin{:});
end