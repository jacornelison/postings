function varargout=arrayfunc(varargin)
% function varargout=arrayfunc(varargin)
%
% Wrapper around arrayfun() which appends "'uniformoutput',false" as arguments
% to arrayfun().
%
  [varargout{1:nargout}] = arrayfun(varargin{:}, 'uniformoutput', false);
end
