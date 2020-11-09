function varargout = paramparse (N, D, varargin)
% [A B C ...] = PARAMPARSE (N, D, args)
%
% Parse lists of parameter-value pairs
%
% N is list of parameter names
% D is list of default values
% args ... are parameter-value pairs
%
% RWO 081211

V = D;
for (ii = 1:2:length(varargin))
	jj = strmatch (varargin{ii}, N, 'exact');
	if length(jj) > 1
		jj = jj(1);
	end;
	if isempty (jj)
		jj = strmatch (varargin{ii}, N);
		if length (jj) > 1
			error (['Parameter name ' varargin{ii} ' has multiple matches.']);
		end;
		if isempty (jj)
			error (['Unknown parameter name ' varargin{ii} '.']);
		end;
	end;
	if length(varargin) > ii
		V{jj} = varargin{ii+1};
	else
		error (['No value given for parameter ' varargin{ii} '.']);
	end;
end;

varargout = V;
