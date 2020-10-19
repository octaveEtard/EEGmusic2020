function x_txt = formatValueDisplay(x,n)
%
% pltools.formatValueDisplay
% Part of the pltools package (github.com/octaveEtard/pltools).
% Author: Octave Etard, 2020
%
if nargin < 2
    n = 1;
end
log_x = floor(log10(abs(x)));
if log_x < -1 || 1 < log_x % use scientific notation
    x_txt = sprintf(['%.',int2str(n),'f\\cdot10^{%i}'],x*10^(-log_x),log_x);
else % ... except for abs(x) from 0.1 to 100
    x_txt = sprintf(['%.',int2str(n),'f'],x);
end
end
%
%