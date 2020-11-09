function tags=list_tags(t1,t2, varargin)
% list_tags(t1,t2, varargin)
% t1 - start time
% t2 - stop time
% varargin: Valid selections include:
%     dk: select a particular dk angle
%     phase: select a particular phase
%     scanset: select a particular scanset
%     has_tod: set to 1 to select only tags w/ tod's
%
% ex: tags=list_tags('2012-Apr-01','2012-Apr-07','dk',113,'phase','B','has_tod');
% rwa - 2012 mar 1

% Default arguments.
phases=[];
dks=[];
scansets=[];
has_tod = 0;

% Parse varargin.
for i=1:length(varargin)
  % phase
  if strcmpi(varargin{i}, 'phase')
    if ischar(varargin{i+1})
      phases = varargin{i+1};
      i = i + 1;
    end
  end

  % dk
  if strcmpi(varargin{i}, 'dk')
    if isnumeric(varargin{i+1})
      dks=varargin{i+1};
      i = i + 1;
    end
  end

  % scanset
  if strcmpi(varargin{i}, 'scanset')
    if isnumeric(varargin{i+1})
      scansets = varargin{i+1};
      i = i + 1;
    end
  end
  
  % has_tod
  if strcmpi(varargin{i}, 'has_tod')
    if isnumeric(varargin{i+1})
      has_tod=varargin{i+1};
    end
  end
end

if has_tod
  tags=get_tags('all','has_tod');
else
  tags=get_tags('all');
end

[days,phase,scanset,dk,extra]=parse_tag(tags);

t1=datenum(t1,'yyyy-mmm-dd');
t2=datenum(t2,'yyyy-mmm-dd');

days=num2str(days');

for ii=1:length(tags)
  dayn(ii)=datenum(days(ii,:),'yyyymmdd');
end

c1=(dayn >= t1 & dayn <= t2);

all=dayn==dayn;

if ~isempty(dks)
  c2=dk==dks;
else
  c2=all;
end

if ~isempty(scansets)
  c3=scanset==scansets;
else
  c3=all;
end

if ~isempty(phases)
  for ii=1:length(tags)
    c4(ii)=strcmpi(phase(ii),phases);
  end
else
  c4=all;
end

cc=c1 & c2 & c3 & c4;
tags=tags(cc);

disp(['list_tags: You have downselected to ' num2str(length(tags)) ' tags']);
