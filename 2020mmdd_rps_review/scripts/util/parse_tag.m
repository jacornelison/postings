function [day,phase,scanset,dk,extra]=parse_tag(tag)
%[day,phase,scanset,dk,extra]=parse_tag(tag)
%
% if tag is '20110101Db01_dk293' day is numeric 20060101
% phase is Db, scanset is 01, deck angle is 293
%
% can handle scanset or phase tags
%
% e.g. for tag='20110521Db01_dk293'
% day=20110521
% phase='Db'
% scanset=1
% dk=293
%

if ~iscell(tag)
  tag={tag};
end

% Grab everything using a complicated regular
% expression
pat=['^' ...                              % start of str
     '(?<tagdate>[0-9]+)' ...             % numeric date
     '(?<tagphase>[A-Z][abcdefg]?)' ...   % 1 or 2 char phase
     '(?<tagsset>[x0-9][0-9])?' ...       % scanset number
     '(_dk)?' ...                         % look for "_dk"
     '(?<tagdk>(?(4)[0-9]+))' ...         % if found, get dk
     '(_)?' ...                           % look for extra "_"
     '(?<tagextra>(?(6).*))$'];           % if found, get extra

% and parse them all at once.
tok=regexp(tag,pat,'names');

N = numel(tag);
day     = NaN(1, N);
phase   = cell(1, N);
scanset = NaN(1, N);
dk      = NaN(1, N);
extra   = cell(1, N);

for i=1:N
  if isempty(tok{i})
    error(['Could not parse tag ' tag{i} '.']);
  end

  day(i)=str2num(tok{i}.tagdate);
  phase{i}=tok{i}.tagphase;
  if isempty(tok{i}.tagsset)
    scanset(i)=NaN;
  elseif tok{i}.tagsset(1)=='x'
    scanset(i)=-1*str2num(tok{i}.tagsset(2:end));
  else
    scanset(i)=str2num(tok{i}.tagsset);
  end
  if isempty(tok{i}.tagdk)
    dk(i)=NaN;
  else
    dk(i)=str2num(tok{i}.tagdk);
  end
  extra{i}=tok{i}.tagextra;
end

if N == 1
  phase = phase{1};
  extra = extra{1};
end

return
