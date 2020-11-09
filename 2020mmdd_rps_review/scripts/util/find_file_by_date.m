function data=find_file_by_date(reqdate,filebase,usecache)
%data=find_file_by_date(reqdate,filebase,usecache)
%
% Read a parameters CSV file with base name filebase whose appended datestamp
% is after reqdate but before the datestamp of the next available file.
%
% INPUTS
%   reqdate     A YYYYMMDD-form integer date. The input file version is chosen
%               to be after this date and before the next available date.
%
%   filebase    A path/file base name. The wildcard '_*.csv' is appended to
%               list all date stamped versions of the file.
%
%   usecache    Defaults to false. If true, the directory listings and any
%               subsequent read data are cached in persistent variable across
%               invocations. Setting to false at any point clears the cache.
%
% OUTPUTS
%   data        A struct of data returned in ParameterRead format.
%

persistent persist_d;
persistent persist_data;

if ~exist('usecache','var') || isempty(usecache)
  usecache = true;
end

% Initialize the caches if they don't exist or clear them out if caching is
% not requested (essentially, allows a reset option)
if ~usecache || isempty(persist_d)
  persist_d = struct('key',{{}}, 'val',{[]});
end
if ~usecache || isempty(persist_data)
  persist_data = struct('key',{{}}, 'val',{[]});
end

% check for a cached directory listing
key = find(strcmp(filebase, persist_d.key));
if ~isempty(key)
  d = persist_d.val{key};
else
  % otherwise get the directory listing and cache that new info
  d = dir(sprintf('%s_*.csv',filebase));

  % Only store the data if a cache is desired
  if usecache
    persist_d.key{end+1} = filebase;
    persist_d.val{end+1} = d;
  end
end
if isempty(d)
  warning('No files matching pattern %s_*.csv', filebase);
  data = [];
  return
end

[pathstr fname]=fileparts(filebase);
for k=1:length(d)
  ind(k)=~isempty(regexp(d(k).name,[fname '_\d{8,8}.csv']));
end
dl=d(ind);

for i=1:length(dl)
  dat_date(i)=str2num(dl(i).name(end-11:end-4));
end

% find the last file whose date is less than the specified date
dlt=find(dat_date<=reqdate);
dlt=dlt(end);

filename=sprintf('%s_%08d.csv',filebase,dat_date(dlt));

% again check for a cached result
key = find(strcmp(filename, persist_data.key));
if ~isempty(key)
  data = persist_data.val{key};
else
  disp(sprintf('find_file_by_date: Reading file %s',filename));
  data=ParameterRead(filename);

  if usecache
    persist_data.key{end+1} = filename;
    persist_data.val{end+1} = data;
  end
end

