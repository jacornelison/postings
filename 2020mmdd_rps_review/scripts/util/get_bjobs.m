% b=get_bjobs
%
% Run bjobs and parse the input into a form
% Matlab can easily use.  Any arguments are
% passed through to bjobs.
%
% Warning: some arguments to bjobs can mess
% up the column alignment and will not work
% with get_bjobs.  These are: -l, -w and -W.
% Note that the option -WL apparently works.
function b=get_bjobs(varargin)

% accept arbitrary arguments, as a string or
% a cell array of strings.  We don't care what
% they are: just feed them to bjobs.
if nargin==0
  args={};
else
  args=varargin;
end
if length(args)==1 && iscell(args{1})
  args=args{1};
end
argstr='';
for i=1:length(args)
  argstr=[argstr ' ' args{i}];
end

% call it!
[s,r]=system(['bjobs' argstr]);

% get first (header) row
headline=strtok(r,10);
% if might have an error message if there's no such job
if ~isempty(strfind(headline,'is not found'))
  b=[];
  return
end

% scan it for column names and sizes
tmp=regexp(headline,'\w+','start');
colidx=[tmp(:) [tmp(2:end)'-1; length(headline)]];
scanfmt='';
for i=1:size(colidx,1)
  colname{i}=lower(strtrim(headline(colidx(i,1):colidx(i,2))));
  if i==size(colidx,1)
    scanfmt=[scanfmt '%[^\n]\n'];
  else
    scanfmt=[scanfmt '%' num2str(colidx(i,2)-colidx(i,1)+1) '[^\n]'];
  end
end

% parse succeeding lines accordingly
P=textscan(r,scanfmt,'HeaderLines',1,'Delimiter','','WhiteSpace','');

% put everything together into output structure
b=[];
for i=1:length(colname)
  b.(colname{i})=P{i};
end

return
