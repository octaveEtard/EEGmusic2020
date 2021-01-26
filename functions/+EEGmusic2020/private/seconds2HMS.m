function t = seconds2HMS(timeInSeconds, formatOut)
% seconds2HMS Convert input 'timeInSeconds' in hours / minutes /
% seconds, and optionally returns it as a string.
%
% Usage: t = seconds2HMS(timeInSeconds, ...)
%
% Required input:
%   timeInSeconds   time (in second) to convert
%
% Optional inputs:
%   formatOut how to format the output string
%
% If formatOut is not set, an array of size 3 with the first element
% containing the hours, the second the minutes and the third the seconds is
% returned, e.g. seconds2HMS(3666) returns [1;1;6].
%
% If formatOut is a char, it will be inserted between h/mn/s, e.g.:
% seconds2HMS(3666, ':') returns '01:01:06'
%
% If formatOut is a cell of strings, it must be of size 3, and the first
% string will be inserted after the hours, the second after the minutes and
% the third after the seconds, e.g.: 
% seconds2HMS(3666, {'h','mn','s'}) returns '01h01mn06s'
%
%

if isnan(timeInSeconds) || timeInSeconds < 0
    t = nan;
    return;
end

% Conversion to hours / mn / s format
t = fix(mod(timeInSeconds, [0; 3600; 60]) ./ [3600; 60; 1]);

% formating of the ouput
if nargin > 1
    t = sprintf('%02i%02i%02i', t);
    if ischar(formatOut)
        t = [t(1:2), formatOut, t(3:4), formatOut, t(5:6)];
    elseif iscell(formatOut) && length(formatOut) == 3
        t = [t(1:2), formatOut{1}, t(3:4), formatOut{2}, t(5:6), formatOut{3}];
    else
        error('"formatOut" input not understood');
    end
end
end