function l = addComparisonLine(ax,x,y,dy,varargin)
%
% pltools.shadeYarea
% Part of the pltools package (github.com/octaveEtard/pltools).
% Author: Octave Etard, 2020
%
if numel(dy) == 1
    dy = [dy,dy];
end

x = [x(1),x,x(2)];
y = [y-dy(1),y,y,y-dy(2)];

l = plot(ax,x,y,varargin{:});

end