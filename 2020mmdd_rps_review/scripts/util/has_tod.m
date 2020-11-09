% function has_tod(tags,datadir,filetype)
% Check whether TOD files exist for the specified tags.
% By default: datadir is 'data/real'
%             filetype is '_tod.mat'
%
% The whole point of this function is to use /bin/ls
% so we don't get hung up forever doing slow directory
% listings on NFS.
function s=has_tod(tags,datadir,filetype)

if ~iscell(tags)
  tags={tags};
end
if nargin<2 || isempty(datadir)
  datadir='data/real/';
end
if nargin<3 || isempty(filetype)
  filetype='_tod.mat';
end

% TOD files are now organized by monthly directories, YYYYMM.
% Construct the month name for each tag.
tagchar=char(tags);
tagmonth=tagchar(:,1:6);
[tagmonth, i, j]=unique(cellstr(tagmonth));

% File name for each tag.
% fnchar=strcat(tagchar,filetype);
fncell=strcat(tags,filetype);
s=false(length(j),1);

% Loop over months.
for i=1:length(tagmonth)

  % Get full listing using /bin/ls.  No wildcards!
  [stat flist]=system_safe(['/bin/ls ' fullfile(datadir,tagmonth{i},'/')]);
  if stat~=0
    error(['Data directory ' fullfile(datadir,tagmonth{i}) ' cannot be read.']);
  end
  if isempty(flist)
    continue
  end

  % And check which tags' files are missing from the list.
  flist=textscan(flist,'%s');
  flist=flist{1};
  % [missing,k]=setdiff(fnchar(j==i,:),flist);
  [missing,k]=setdiff(fncell(j==i),flist);
  tmps=true(sum(j==i),1);
  tmps(k)=false;

  % Put back in to overall mask.
  s(j==i)=tmps;

end

return
