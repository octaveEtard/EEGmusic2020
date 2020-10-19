function varargout = putInCell(varargin)
%
% Put input arguments in cell if that's not already the case
%
for iarg = 1:nargin
    if ~iscell(varargin{iarg})
        varargin{iarg} = varargin(iarg);
    end
end

varargout = varargin;

end
%
%