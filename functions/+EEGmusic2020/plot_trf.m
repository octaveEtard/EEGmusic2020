function [sgHandles,nsHandles] = plot_trf(ax,x,y,err,ROI,sg,...
    sgArgs,nsArgs,shadingArgs,...
    col,col_null,collShading,...
    lsty,fts,lgValues)
%
% EEGmusic2020.plot_trf
% Part of the EEGmusic2020 code.
% Author: Octave Etard
%
% Plot mean TRFs curve and shading around the mean.
%
n = size(y,2);

if size(sg,2) == 1
    sg = repmat(sg,[1,n]);
end

% for the legend
nsHandles = gobjects(n,1);
sgHandles = gobjects(n,1);

% for each column in y, plot y(:,iCol) against x, with a shading
% corresponding to +/- err(:,iCol) around y(:,iCol)
for itc = 1:n
    shadingArgs_ = [shadingArgs,'FaceColor',collShading(itc,:)];
    nsArgs_ = [nsArgs,'Color',col_null(itc,:),'LineStyle',lsty{itc}];
    sgArgs_ = [sgArgs,'Color',col(itc,:),'LineStyle',lsty{itc}];
    
    [nsHandles(itc),sgHandles(itc)] = pltools.plot_tc(ax,x,y(:,itc),...
        err(:,itc),ROI,sg(:,itc),nsArgs_,sgArgs_,shadingArgs_);
end

% make sure trf curves are on top of shading
if ~isa(nsHandles,'matlab.graphics.GraphicsPlaceholder')
    uistack(nsHandles,'top');
end
if isa(sgHandles,'matlab.graphics.GraphicsPlaceholder')
    lgHandle = nsHandles;
else
    uistack(sgHandles,'top');
    lgHandle = sgHandles;
end

ax.XAxis.Label.String = 'Time (ms)';
ax.YAxis.Label.String = 'Amplitude (a.u.)';

if ~isempty(lgValues)
    legend(lgHandle,lgValues,...
        'box','off',...
        'Location','northeast',...
        'FontSize',fts);
end
m = max(abs([y+err,y-err]),[],'all');
ax.YAxis.Limits = m*[-1,1];
end
%
%