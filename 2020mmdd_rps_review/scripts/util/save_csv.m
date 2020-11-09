% save_csv(fname,p) where p is a structure
% saves the fields of p in a CSV file.  This
% is a wrapper around ParameterWrite that
% primarily makes it easier to generate the
% necessary header with units and data types.
%
% save_csv(...,'comment','Example comment')
% fills in the specified comments.  Can be
% a string or cell array of strings.
%
% save_csv(...,'units',u) specifies the units
% for the fields, where u should be a
% structure with the same fields as p.
%
% Example:
%  p = [];
%  p.gcp = (1:2640)-1;
%  p.opteff = opteff;
%  save_csv('opteff.csv',p)
function save_csv(fname,p,varargin)

cmt = '';
u = [];

for i=1:2:length(varargin)
  s = varargin{i};
  if ~ischar(s) || isempty(s)
    error(['Invalid data type for parameter name.']);
  end
  if length(varargin)<(i+1)
    error(['No value specified for parameter ' s '.']);
  end
  val = varargin{i+1};
  switch(lower(s))
    case {'comment','cmt'}, cmt = val;
    case {'unit','units'}, u = val;
    otherwise, error(['Unknown parameter ' s '.']);
  end
end

if isempty(cmt)
  cmt = ['# ' datestr(now)];
end
if ischar(cmt)
  cmt = {cmt};
end
for i=1:length(cmt)
  if isempty(cmt{i}) || cmt{i}(1)~='#'
    cmt{i} = ['# ' cmt{i}];
  end
end

k = [];

k.comments = cmt;
k.fields = fieldnames(p);
for i=1:length(k.fields)
  v = p.(k.fields{i});
  if iscell(v)
    if ischar(v{1})
      k.formats{i} = 'string';
    elseif isnumeric(v{1})
      k.formats{i} = 'double';
    end
  elseif isinteger(v) || all(v - floor(v) == 0)
    k.formats{i} = 'integer';
  else
    k.formats{i} = 'double';
  end
  if isstruct(u) && isfield(u,k.fields{i})
    k.units{i} = u.(k.fields{i});
  elseif strcmp(k.formats{i},'string')
    k.units{i} = '-';
  else
    k.units{i} = '(null)';
  end
  if ~isvector(v)
    error(['Field ' k.fields{i} ' is not a vector.']);
  end
end

ParameterWrite(fname,p,k);

return
