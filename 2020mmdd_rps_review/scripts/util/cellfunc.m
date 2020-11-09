function varargout=cellfunc(varargin)
% function varargout=cellfunc(varargin)
%
% Wrapper around cellfun() which appends "'uniformoutput',false" as arguments
% to cellfun().
%
  [varargout{1:nargout}] = cellfun(varargin{:}, 'uniformoutput', false);
end
