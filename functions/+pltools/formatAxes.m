function axes = formatAxes(axes, xOffset, yOffset, width, height)
%
% pltools.formatAxes
% Part of the pltools package (github.com/octaveEtard/pltools).
% Author: Octave Etard, 2020
%
% formatAxes Format axes in 'axes' so that the tight bounding box of the
% ensemble fits the rectangle with bottom left corner (xOffset, yOffet),
% and dimension 'width' x 'height' (normalised values)

nAxes = length(axes);

for iAx = 1:nAxes
    if isprop(axes(iAx),'ActivePositionProperty')
        axes(iAx).ActivePositionProperty = 'position';
    end
    axes(iAx).Units = 'normalized';
end


%% left, bottom, right, top coordinate of all axes
left = zeros(nAxes, 1);
bottom = zeros(nAxes, 1);
right = zeros(nAxes, 1);
top = zeros(nAxes, 1);

for iAx = 1:nAxes
    left(iAx) = axes(iAx).Position(1);
    bottom(iAx) = axes(iAx).Position(2);
    right(iAx) = axes(iAx).Position(1) + axes(iAx).Position(3);
    top(iAx) = axes(iAx).Position(2) + axes(iAx).Position(4);
end


%% determining the bounding box of the plot and the distances to scale (axis)
% distance made of axes on the left of each axis -- will be scaled
hScalable = zeros(nAxes, 1);
% distance made of axes below  each axis -- will be scaled
vScalable = zeros(nAxes, 1);

[~, isLeft] = sort(left);
rightMost = right(isLeft(1));
hSize = rightMost - left(isLeft(1));

[~, isBottom] = sort(bottom);
topMost = top(isBottom(1));
vSize = topMost - bottom(isBottom(1));

% bounding box coordinates
[bbLeft,bbBottom,bbRight,bbTop] = getBBox(axes(1));

% bbLeft = axes(1).Position(1) - axes(1).TightInset(1);
% bbBottom = axes(1).Position(2) - axes(1).TightInset(2);
% bbRight = axes(1).Position(1) + axes(1).Position(3) + axes(1).TightInset(3);
% bbTop = axes(1).Position(2) + axes(1).Position(4) + axes(1).TightInset(4);

for iA = 2:nAxes
    r = right(isLeft(iA));
    l = left(isLeft(iA));
    hScalable(isLeft(iA)) = hSize - max(0, rightMost - l);
    
    if r > rightMost
        hSize = hSize + r - l - max(rightMost - l, 0);
        rightMost = r;
    end
    
    t = top(isBottom(iA));
    b = bottom(isBottom(iA));
    vScalable(isBottom(iA)) = vSize - max(0, topMost - b);
    if t > topMost
        vSize = vSize + t - b - max(topMost - b, 0);
        topMost = t;
    end
    
    [bbLeftTmp,bbBottomTmp,bbRightTmp,bbTopTmp] = getBBox(axes(iA));
    bbLeft = min(bbLeft, bbLeftTmp);
    bbBottom = min(bbBottom, bbBottomTmp);
    bbRight = max(bbRight, bbRightTmp);
    bbTop = max(bbTop, bbTopTmp);
    
    % ax = axes(iA);
    % bbLeft = min(bbLeft, ax.Position(1) - axes(iA).TightInset(1));
    % bbBottom = min(bbBottom, ax.Position(2) - axes(iA).TightInset(2));
    % bbRight = max(bbRight, ax.Position(1) + ax.Position(3) + ax.TightInset(3));
    % bbTop = max(bbTop, ax.Position(2) + ax.Position(4) + ax.TightInset(4));
end

% % scaling factors
% w = (1 - (bbRight - bbLeft +(1-width)) + hSize) / hSize
% h = (1 - (bbTop - bbBottom +(1-height)) + vSize) / vSize
%
% for iAx = 1:nAxes
%     p = axes(iAx).Position;
%     % translation
%     p = [p(1) - bbLeft + xOffset, p(2) - bbBottom + yOffset, p(3), p(4)];
%     % scaling
%     p = [p(1)+hScalable(iAx)*(w-1), p(2)+vScalable(iAx)*(h-1), p(3)*w, p(4)*h];
%     % offset
%     axes(iAx).Position = [p(1), p(2), p(3), p(4)];
% end

% scaling factors
w = (width - (bbRight - bbLeft) + hSize) / hSize;
h = (height - (bbTop - bbBottom) + vSize) / vSize;

for iAx = 1:nAxes
    p = axes(iAx).Position;
    % translation
    p = [p(1) - bbLeft + xOffset, p(2) - bbBottom + yOffset, p(3), p(4)];
    % scaling
    p = [p(1)+hScalable(iAx)*(w-1), p(2)+vScalable(iAx)*(h-1), p(3)*w, p(4)*h];
    % offset
    axes(iAx).Position = [p(1), p(2), p(3), p(4)];
    moveBackText(axes(iAx),w,h);
end



end
%%
function [bbLeft,bbBottom,bbRight,bbTop] = getBBox(ax)

ax.Units = 'normalized';

xAx = ax.Position(1);
yAx = ax.Position(2);
widthAx = ax.Position(3);
heightAx = ax.Position(4);

txt = findobj(ax.Children,'Type','Text');

% ax bounding box
switch ax.Type
    case 'axes'
        bbLeft = xAx - ax.TightInset(1);
        bbBottom = yAx - ax.TightInset(2);
        bbRight = xAx + widthAx + ax.TightInset(3);
        bbTop = yAx + heightAx + ax.TightInset(4);
    case 'colorbar'
        bbLeft = xAx;
        bbBottom = yAx;
        bbRight = xAx + widthAx;
        bbTop = yAx + heightAx;
        
        % Hack: using colorbar label to find tick extent
        txt = [text,ax.Label];
end


nTxt = length(txt);

for iTxt = 1:nTxt
    txt(iTxt).Units = 'normalized';
    % in the parent axis coordinates
    x = txt(iTxt).Extent(1);
    y = txt(iTxt).Extent(2);
    width = txt(iTxt).Extent(3);
    height = txt(iTxt).Extent(4);
    % in the figure coordinate
    x = xAx + x*widthAx;
    y = yAx + y*heightAx;
    width = width * widthAx;
    height = height * heightAx;
    
    bbLeft = min(bbLeft,x);
    bbRight = max(bbRight,x+width);
    
    bbTop = max(bbTop,y+height);
    bbBottom = min(bbBottom,y);
    
end
end

function moveBackText(ax,w,h)

txt = findobj(ax.Children,'Type','Text');
nTxt = length(txt);

for iTxt = 1:nTxt
    txt(iTxt).Units = 'normalized';
    % in the parent axis coordinates
    x = txt(iTxt).Position(1);
    y = txt(iTxt).Position(2);
    
    if x < 0.5
        refX = 0;
    else
        refX = 1;
    end
    if y < 0.5
        refY = 0;
    else
        refY = 1;
    end
    
    txt(iTxt).Position(1) = refX + (x-refX)/w;
    txt(iTxt).Position(2) = refY + (y-refY)/h;
end


end
%
%