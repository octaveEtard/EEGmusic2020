function formatAxisLabels(graphObj,varargin)
%
% pltools.formatAxisLabels
% Part of the pltools package (github.com/octaveEtard/pltools).
% Author: Octave Etard, 2020
%
%% default values
ip = inputParser;

ip.addRequired('graphObj', @isGraphicsArray);

ip.addOptional('ftsLabel',11,@isnumeric);
ip.addOptional('lwd',0,@isnumeric);
ip.addOptional('interpreter','tex',@ischar);
ip.addOptional('ftsTick',0,@isnumeric);
ip.addOptional('ftsTitle',0,@isnumeric);


ip.parse(graphObj,varargin{:});

ftsLabel = ip.Results.ftsLabel;
lwd = ip.Results.lwd;
interpreter = ip.Results.interpreter;
ftsTick = ip.Results.ftsTick;
ftsTitle = ip.Results.ftsTitle;


if ftsTitle == 0
    ftsTitle = ftsLabel;
end
if ftsTick == 0
    ftsTick = round(0.8 * ftsLabel);
end

set(0,'defaultTextInterpreter',interpreter);

%%
n = numel(graphObj);
for iAx = 1:n
    ax = graphObj(iAx);
    ax.Box = 'off';
    
    % set all to bet displayed with latex
    
    ax.TickLabelInterpreter = interpreter;
    ax.XLabel.Interpreter = interpreter;
    ax.YLabel.Interpreter = interpreter;
    ax.Title.Interpreter = interpreter;
    
    try
        % this call will also affect labels, title, legend...
        % hence it should come first, not to override the later calls.
        ax.FontSize = ftsTick;
        
        
        ax.XAxis.Label.FontSize = ftsLabel;
        ax.YAxis.Label.FontSize = ftsLabel;
        
        if lwd ~= 0
            ax.XAxis.LineWidth = lwd;
            ax.YAxis.LineWidth = lwd;
        end
    end
    
    ax.Title.FontSize = ftsTitle;
    ax.Title.FontWeight = 'normal';
    
    if isfield(ax,'Legend')
        ax.Legend.FontSize = ftsLabel;
        ax.Legend.Title.FontSize = ftsLabel;
        ax.Legend.Interpreter = interpreter;
    end
      
end

end

function b = isGraphicsArray(x)
b = all(isgraphics(x(:)));
end