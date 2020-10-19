function printFigure(figHandle, pathToFolder, fileName, resolution, widthPaper, heightPaper, w, h, doPrint, format, doScale)
%
% pltools.printFigure
% Part of the pltools package (github.com/octaveEtard/pltools).
% Author: Octave Etard, 2020
%
if nargin< 9
    doPrint = true;
end
if nargin < 10
    format = 'pdf'; %'svg'
end
if nargin < 11
    doScale = true;
end

switch format
    case 'emf'
        format = '-dmeta';
    otherwise
        format = ['-d',format];
end

if isempty(pathToFolder)
    pathToFolder = '.';
end

figHandle.Units = 'centimeters';
figPos = figHandle.Position;
figHandle.Position = [figPos(1), figPos(2), widthPaper, heightPaper];

figHandle.PaperUnits = 'centimeters';
figHandle.PaperSize = [widthPaper, heightPaper];
figHandle.PaperPosition = [0, 0, widthPaper heightPaper];

if doScale
    axes = figHandle.Children(arrayfun(@(x) any(strcmp(x.Type, {'axes','colorbar'})), figHandle.Children));
    pltools.formatAxes(axes, (1-w)/2, (1-h)/2, w, h);
end

% if doPrint
%     figHandle.Color = 'none';
%     arrayfun(@(ax) set(ax,'Color','none'),figHandle.Children);
% end

% % formatAxes(axes, 0, (1-h)/2, w, h);

%% Printing
if doPrint
    [~,fileName,~] = fileparts(fileName);
    pathToFile = fullfile(pathToFolder, fileName);
    res = ['-r' num2str(resolution)];
    
    if ~exist(pathToFolder,'dir')
        mkdir(pathToFolder);
    end
    
    print(figHandle, pathToFile, format, res);
end
end
%
%