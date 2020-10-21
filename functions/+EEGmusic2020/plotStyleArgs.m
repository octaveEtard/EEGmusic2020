function [col,col_null,shadingArgs,nsArgs,sgArgs] = plotStyleArgs()
%
% EEGmusic2020.plotStyleArgs
% Part of the EEGmusic2020 code.
% Author: Octave Etard, 2020
%
% Return some plotting arguments to ensure consistent style across figures.
%

% default colors
col = pltools.getColorOrder(true,3); % (black & red)
col(2,:) = [];
col = [col ; .85, .85, .85]; % and a light grey for lines linking subjects

% lighter versions of above black & red for null distributions
col_null = [...
    153,153,153;...
    230,157,108] / 255;

% shading
shadingArgs = {...
    'FaceAlpha',0.3,...
    'EdgeColor','none'};

% default / ns line thickness
nsArgs = {...
    'Linewidth',1};

% bold for significant regions
sgArgs = {...
    'Linewidth',2};
end
%
%