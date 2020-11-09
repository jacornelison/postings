% ind=has_pairmap(tags,datadir,filetype)
%  or
% ind=has_pairmap(tags,mapopt)
%
% Check whether any pairmap files exist for the specified tags.
% By default: datadir is 'pairmaps/0000/real/'
%             filetype is '_filtp3_weight3_gs.mat'
%
% The whole point of this function is to use /bin/ls
% so we don't get hung up forever doing slow directory
% listings on NFS.
%
% In order to make it easier to use this for setting
% up mapmaking jobs, you can also call as
%   has_pairmap(tags,mapopt)
% in which case the datadir and filetype will be
% automatically guessed for you.
function s=has_pairmap(tags,datadir,filetype)

if ~iscell(tags)
  tags={tags};
end

% see if we were passed a mapopt structure
if nargin==2 && isstruct(datadir)
  mapopt=datadir;
  
  % fill in default mapopts if not specified
  mapopt=get_default_mapopt(mapopt);
  
  sernum=mapopt.sernum;
  if isnumeric(sernum)
    sernum=num2str(sernum,'%.8d');
  end
  datadir=fullfile('pairmaps',sernum(1:4),sernum(5:8));
  weight=mapopt.weight;
  if isnumeric(mapopt.weight)
    weight=num2str(mapopt.weight);
  end
  gsstr='';
  if isfield(mapopt,'gs') && ~isempty(mapopt.gs)
    if ischar(mapopt.gs) && strcmp(mapopt.gs,'1')
      gsstr='_gs';
    elseif mapopt.gs
      gsstr='_gs';
    end
  end
  dpstr='';
  if isfield(mapopt,'deproj') && ~isempty(mapopt.deproj)
    if ischar(mapopt.deproj) && strcmp(mapopt.deproj,'1')
      dpstr='_dp';
    elseif any(mapopt.deproj)
      dpstr='_dp';
    end
  end
  projstr='';
  if isfield(mapopt,'proj') && ~isempty(mapopt.proj)
    if ~strcmp(mapopt.proj,'radec')
      projstr=['_' mapopt.proj];
    end
  end
  filetype=['_filt' mapopt.filt '_weight' weight gsstr dpstr projstr '.mat'];
else
  % datadir, filetype given separately
  if nargin<2 || isempty(datadir)
    datadir='pairmaps/0000/real/';
  end
  if nargin<3 || isempty(filetype)
    filetype='_filtp3_weight3_gs.mat';
  end
end

% File name for each tag.
% fnchar=strcat(tagchar,filetype);
fncell=strcat(tags,filetype);
s=true(length(tags),1);

if ~exist(datadir,'dir')
  s=false(size(s));
  return
end

% Get full listing using /bin/ls.  No wildcards!
[stat flist]=system_safe(['/bin/ls ' fullfile(datadir,'/')]);
if stat~=0
  error(['Pairmap directory ' fullfile(datadir) ' cannot be read.']);
end
if isempty(flist)
  s=false(size(s));
  return
end

% And check which tags' files are missing from the list.
flist=textscan(flist,'%s');
flist=flist{1};
% [missing,k]=setdiff(fnchar(j==i,:),flist);
[missing,k]=setdiff(fncell,flist);

% Put back in to overall mask.
s(k)=false;

return
