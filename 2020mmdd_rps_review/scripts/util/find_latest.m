function latest=find_latest(reqdate, base, ext, direc)
% latest=find_latest(reqdate, base, ext, direc)
%
% Find the date-suffixed file which is the immediately prior or next file
% when compared to the reference date.
%
% INPUTS
%   reqdate    Defaults to datestr(now(), 'yyyymmdd'). Date string as
%              'YYYYMMDD' or decimal integer YYYYMMDD representing the
%              bounding date.
%
%   base       Base name (optional path plus file name prefix) for the files
%              to match.
%
%   ext        Defaults to '.csv'; the file extension to ignore from the
%              suffix.
%
%   direc      Defaults to 'before'. If 'before', gives the latest file with
%              a date suffix not greater than the reference date. If 'after',
%              gives the first file following the reference date.
%
% RETURNS
%   latest     Path to the latest file consistent with the base naming.
%

  if ~exist('reqdate','var') || isempty(reqdate)
    reqdate = datestr(now,'yyyymmdd');
  end
  if ~exist('base','var') || isempty(base)
    base = '';
  end
  if ~exist('ext','var') || isempty(ext)
    ext = '.csv';
  end
  if ~exist('direc','var') || isempty(direc)
    direc = 'before';
  end

  if ischar(reqdate)
    reqdate = str2num(reqdate);
  end
  regext = strrep(ext, '.', '\.');

  % List all matching files
  candids = dir(sprintf('%s*%s', base, ext));
  [path fname] = fileparts(base);
  % Filter the candidate list for exact matches to the form
  %     {base}YYYYMMDD{ext}
  matches = regexp({candids(:).name}, ['^' fname '\d{8,8}' regext '$']);
  matches = cellfun(@(l) ~isempty(l), matches);
  candids = {candids(matches).name};

  % Get the dates for the remaining files
  dates = cellfun(@(c) str2num(c(end-11:end-4)), candids);
  % Guarantee ordering
  [dates,order] = sort(dates);
  candids = candids(order);

  switch lower(direc)
    case 'before'
      candids = candids(dates <= reqdate);
      idx = length(candids);
    case 'after'
      candids = candids(dates > reqdate);
      idx = 1;
    otherwise
      error('unrecognized direction ''%s''', direc)
  end

  if isempty(candids)
    latest = [];
  else
    latest = fullfile(path, candids{idx});
  end

end

