function [ns,sg,shading] = plot_tc(ax,x,y,err,ROI,sg,...
    nsArgs,sgArgs,shadingArgs)
%
% pltools.plot_tc
% Part of the pltools package (github.com/octaveEtard/pltools).
% Author: Octave Etard, 2020
%
% Plot y with respect to x in axes ax, with style sgArgs in ROI & where sg
% and style nsArgs otherwise. If ~isempty(err) add a shading  around y 
% at +/- err at each point x.
%
%
% --- plotting shading at +/- 1 error around y first
if ~isempty(err)
    yMin = y - err;
    yMax = y + err;
    shading = shadeYarea(ax,x,yMin,yMax,shadingArgs{:}); % FIXME
end

% --- plot y
if ~isempty(ROI) && any( sg )
    % --- some regions are significant:
    % distinguishing between significant and ns regions in the plot
    
    % first plotting everything with col_null colours
    ns = plot(ax,x,y,nsArgs{:});
    
    % then plotting only the significant parts
    x_roi = x(ROI);
    
    y_roi = y(ROI);
    y_roi( ~sg ) = nan;
    
    sg = plot(ax,x_roi,y_roi,sgArgs{:});
    
else
    % --- no significance, plotting everything with one style (nsArgs)
    sg = gobjects(1);
    
    ns = plot(ax,x,y,nsArgs{:});
end
end
%
%