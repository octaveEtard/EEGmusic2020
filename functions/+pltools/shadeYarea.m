function p = shadeYarea(ax,x,yMin,yMax,varargin)
%
% pltools.shadeYarea
% Part of the pltools package (github.com/octaveEtard/pltools).
% Author: Octave Etard, 2020
%
% Shade area alongside x, between yMin and yMax
% e.g. to represent error along y for all x.
%
% Can take all optional arguments from patch as options.
%
assert(length(yMin) == length(yMax));

if length(x) ~= length(yMin)
    assert(length(yMin) == 1);
    x = [x(1),x(end)];
    yMin = [yMin,yMin];
    yMax = [yMax,yMax];
end
% row vector
x = x(:)';
yMin = yMin(:)';
yMax = yMax(:)';

x = [x,fliplr(x)];
y = [yMin(:)',fliplr(yMax(:)')];
 
p = patch(ax,'XData',x,'YData',y,varargin{:});
%
end