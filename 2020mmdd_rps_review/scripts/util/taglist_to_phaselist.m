function [phases,idx]=taglist_to_phaselist(tags,varargin)
% phases=taglist_to_phaselist(tags)
%
% For BICEP2 reduction is per 1 hour chunk. This function reduces the
% tag list to a list of tags in each phase
%
% Assumed that tag names have form ...XX_...
% where ... is anything and XX is the scanset number
%
% taglist_to_phaselist(tags,'inclusive') tries to
% join incomplete phases after a restart, as long as
% the dk angles are the same, and the start days are
% not too far apart.
%
% taglist_to_phaselist(tags,'raster') joins together
% phases that are in the same "raster": B+C, D+E+F,
% G+H+I.  You can use 'raster' and 'inclusive'
% together.
%
% [ph,idx]=taglist_to_phaselist(tags) also returns an
% index idx such that tags{i} is in phase idx(i).

if(isempty(tags))
  phases.name=[];
  phases.tags=[];
  return
end

if(~iscell(tags));
  tags=cellstr(tags);
end

% By default, don't join any phases together
do_inclusive=false;
do_raster=false;

% Parse any extra option strings
if nargin>1
  for i=1:length(varargin)
    if ~ischar(varargin{i})
      error(['Don''t know what to do with non-string option.']);
    end
    switch(lower(varargin{i}))
      case 'inclusive', do_inclusive=true;
      case 'raster', do_raster=true;
      otherwise, error(['Unrecognized option "' varargin{i} '".']);
    end
  end
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
tok=regexp(tags,pat,'names');

% validate.
for i=1:length(tok)
  if isempty(tok{i})
    error(['Failed to parse tag ' tags{i} '.']);
  end
end
tok=[tok{:}];

% datenum for schedule start days
dn=datenum({tok(:).tagdate},'yyyymmdd');
dn=dn(:);

% To avoid costly str2num operations, use 'unique' to just
% index all the unique values for a given field.  This works
% for phase and dk, not for day (since we need a numeric
% comparison of days).

% phase
[phtmp,phi,phnum]=unique({tok(:).tagphase});
phnum=phnum(:);

% dk
[dktmp,dki,dknum]=unique({tok(:).tagdk});
dknum=dknum(:);

% In different versions of Matlab, dn, phnum, and dknum can each have different
% shapes (i.e. row or vector), so rather than manually transposing them to a
% particular shape, let Matlab choose a consistent shape by using the (:)
% syntax so that the binary-or below will succeed.

% Find tags that have different schedule day, dk, or
% phase number
cnew=false(numel(dn),1);
cnew(2:end)=(diff(dn)~=0 | diff(dknum)~=0 | diff(phnum)~=0);
cnew(1)=true;

% Join incomplete phases if requested, i.e. for deproj
if do_inclusive
  last_dn=dn(1);
  for i=find(cnew)'
    if i==1
      continue
    end
    % Join if compatible.  Could also check scanset(i)>scanset(i-1)
    % but I have left this out for now.  Note that dk can be NaN
    % if some tags don't have a dk angle.
    if (phnum(i)==phnum(i-1) ...
        && dknum(i)==dknum(i-1) ...
        && dn(i)-last_dn<3)
      cnew(i)=false;
    else
      last_dn=dn(i);
    end
  end
end

% Join phases together in rasters, i.e. for deproj 'by_raster'
if do_raster
  for i=find(cnew)'
    if i==1
      continue
    end
    % Join if not a "raster-end" phase.  This is lame... matches
    % current logic in reduc_coaddpairmaps.
    if ~ismember(tok(i).tagphase,'CFI')
      cnew(i)=false;
    end
  end
end

% Having identified first tag of each phase, construct output
p={};
t={};
idx=zeros(size(cnew));
jnew=find(cnew);
for i=1:length(jnew)
  p{i}=[tok(jnew(i)).tagdate tok(jnew(i)).tagphase];
  if ~isempty(tok(jnew(i)).tagdk)
    p{i}=[p{i} '_dk' tok(jnew(i)).tagdk];
  end
  if i==length(jnew)
    t{i}=tags(jnew(i):end);
    idx(jnew(i):end)=i;
  else
    t{i}=tags(jnew(i):(jnew(i+1)-1));
    idx(jnew(i):(jnew(i+1)-1))=i;
  end
end

phases.name=p;
phases.tags=t;

return
