function d = readarc_remote (host, varargin)
% READARC_REMOTE  extracts fields from arc files on a remote
%  machine.  Do not use readard_remote directly.  Instead,
%  use load_arc and put the host (and optional username) in
%  the path, like this:
%
%  d = load_arc ('reuben@bicep0.caltech.edu:/data/bicep2daq/arc/', ...)
%

% Default port - change if needed
port = 22;

% Parse host name, possibly including username
[h1 h2] = strtok (host, '@');
if ~isempty (h2)
	username = h1;
	host = h2(2:end);
else
	username = getenv ('USER');
end;
if isempty (username)
	username = 'bicep';
end;

test_host = host;
if strcmp(test_host,'localhost') || strcmp(test_host,'127.0.0.1')
	test_host = getenv('HOSTNAME');
end

% Set path to the executable appropriately for the host
switch(test_host)
	case {'spud','spud.spa.umn.edu','spud.spa'},
		dump_exec = '/home/reuben/tes_analysis/mex_code/arcfile/dumparc';
	otherwise,
		dump_exec = '/home/bicep0/reuben/tes_analysis/mex_code/arcfile/dumparc';
end;

% Decide where to put temporary file
tmppath = strtok (userpath, pathsep);
if isempty (tmppath) || ~exist (tmppath, 'dir')
	tmppath = pwd;
end;
d = dir (fullfile (tmppath, 'tmp_arcfile_buffer_*'));
if ~isempty (d)
	disp (['Found old tmp buffer files in ' tmppath ' - consider deleting them.']);
	for (ii = 1:length (d))
		disp (['     ' d(ii).name]);
	end;
	disp ([' ']);
end;
tmpfile = fullfile (tmppath, ['tmp_arcfile_buffer_' datestr(now,'yyyymmdd_HHMMSS')]);
disp (['Saving data in temporary file ' tmpfile]);

% Assemble the SSH command
ssh_cmd = ['ssh -p ' num2str(port) ' -l ' username ' ' host ' '];
exec_cmd = [dump_exec ' '];
for (ii = 1:length (varargin))
  exec_cmd = [exec_cmd ' ''' varargin{ii} ''''];
end;
exec_cmd = ['"' exec_cmd '"'];

% disp (['command line: ' ssh_cmd exec_cmd ' > ' tmpfile]);

% And run the SSH command
s = unix ([ssh_cmd exec_cmd ' > ' tmpfile]);
if (s ~= 0)
  if exist (tmpfile, 'file')
    delete (tmpfile);
  end;
  error (['Remote process failed!']);
end;
% disp (['Read ' num2str(length(dat)) ' characters over ssh.']);
tmpinfo = dir (tmpfile);
disp (['Read ' num2str(tmpinfo(1).bytes) ' bytes over ssh.']);


% Parse the file
d = [];
f = fopen (tmpfile, 'rt');
got_start = 0;
while ~feof (f)
  if ~got_start
    ll = fgetl (f);
    if ~strcmp (ll, '###')
      continue;
    end;
  end;

  % Line 1 : register name
  ll = fgetl (f);
  [map ll] = strtok (ll);
  [board ll] = strtok (ll);
  [reg ll] = strtok (ll);

  % Line 2 : data type
  ll = fgetl (f);
  tt = ll;

  % Line 3 : size (samples x channels)
  ll = fgetl (f);
  sz = str2num (ll);

  % Initialize output array
  d.(map).(board).(reg) = zeros (sz(1), sz(2), tt);

  % Following lines : data
  ll = '';
  j = 0;
  n = 0;
  ll = fgetl (f);
  while (~feof (f) && ~strcmp (ll, '###'))
    ll = ll - '0';
    ll (ll > 9) = ll (ll > 9) + '0' - 'a' + 10;
    tmp1 = ll (1:2:end);
    tmp2 = ll (2:2:end);
    if (length (tmp1) > length (tmp2))
      tmp1 = tmp1 (1:length(tmp2));
    end;
    try
      tmpb = typecast (uint8(tmp1 * 16 + tmp2), tt);
      n = length (tmpb);
      d.(map).(board).(reg)(j + (1:n)) = tmpb;
      j = j + n;
    catch, keyboard; end;
    ll = fgetl (f);
    % disp ('OK.');
  end;
  got_start = (strcmp (ll, '###'));
end;
fclose (f);
% If something went wrong, print out any output
if isempty (d)
	type (tmpfile);
end;
disp (['Deleting temporary file ' tmpfile]);
delete (tmpfile);
