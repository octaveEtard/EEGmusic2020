function txt = topLeftLabel(ax,label,fts)
%
% pltools.topLeftLabel
% Part of the pltools package (github.com/octaveEtard/pltools).
% Author: Octave Etard, 2020
%
if nargin < 3 
    fts = ax.Title.FontSize;
end

% place the label in the topleft corner of the tight inset of ax
ax.Units = 'normalized';

% ax position is with respect to figure, but text position is with respect
% to ax!
% converting from fig coordinate to ax coordinate
x = - ax.TightInset(1) / ax.Position(3);
y = 1 + ax.TightInset(4) / ax.Position(4);

txt = text(ax,x,y,label,...
    'Units','normalized',...
    'FontSize',fts,...
    'FontWeight','bold',...
    'VerticalAlignment','cap',...
    'HorizontalAlignment','left');

end