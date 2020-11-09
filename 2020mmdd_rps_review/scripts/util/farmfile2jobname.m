function jobname=farmfile2jobname(farmfile,verify)
% jobname=farmfile2jobname(farmfile,verify)
%
% Transforms a farm file path name into a babysitjobs-compatible job name.
%
% INPUTS
%   farmfile    A farm file path string. One compatible form is the standard
%               form 'farmfiles/NNNN/XXXY_ID' where NNNN is a 4-digit numeric
%               serial number, XXXY is either a simulation number and type or
%               the string literal 'real', and ID is any arbitrary string
%               (usually chosen to identify the operation specifically). The
%               second is 'farmfiles/ID' where ID is again an arbitrary
%               string.
%
%   verify      Defaults to false. If true, the farm file is verified to
%               exist, and if it does not, an error is emitted. Also verifies
%               that the farmfile adhere's to the expected conventions.
%
% OUTPUTS
%   jobname     A job name string. babysitjobs has a corresponding structure
%               which takes the form 'NNNN_XXXY_ID' and 'BICEP_ID'.
%
% SEE ALSO
%   farmfilejobname
%   jobname2farmfile

  if ~exist('verify','var') || isempty(verify)
    verify = false;
  end

  is4num = @(s) strcmp(s, sprintf('%04d',str2num(s)));

  % (At least try) to take apart the farmfile path
  [base,id,ext] = fileparts(farmfile);

  % Detect BICEP_* jobs by recognizing that it doesn't have a serial number
  % layout just like babysitjobs() does.
  [dum,maybesernum] = fileparts(base);
  if ~is4num(maybesernum)
    isbicep = true;
    % Keep any extra path components in base (i.e. those after 'farmfiles/')
    % as a prefix which must be handled.
    slashidx = find(base == filesep(), 1, 'first');
    if ~isempty(slashidx)
      prefix = [base(slashidx+1:end) filesep()];
      base = base(1:slashidx);
    else
      prefix = '';
    end
  else
    isbicep = false;
  end

  % Split off the deepest directory name from the rest of the base path
  % Split again to get the serial number separated from the rest of the base
  [base,sernum] = fileparts(base);
  if ~isbicep && isempty(sernum)
    warning('farmfilejobname:nameFormat', ...
        'No serial number could be identified.')
  end

  if verify
    % Does it use farmfiles/ in the beginning?
    hasfarmfiles = strncmp(base, 'farmfiles', length('farmfiles'));
    % Does the serial number match 4 digit number?
    hassernum = isbicep || is4num(sernum);
    % Does the rest of the serial number follow in the id?
    hassimnum = isbicep || strcmp(id(1:4), 'real') || is4num(id(1:4));

    if ~hasfarmfiles || ~hassernum || ~hassimnum
      error('farmfilejobname:nameFormatError', ...
        'farmfile `%s` does not conform to expected format', farmfile)
    end

    % Check that the file exists last since it's presumably the most
    % expensive operation to wait for the file system.
    if ~exist(farmfile)
      error('farmfilejobname:nameFormatError', ...
          'farmfile `%s` does not exist', farmfile)
    end
  end

  % Put all the pieces back together for a job name
  if isbicep
    jobname = sprintf('%sBICEP_%s', prefix, id);
  else
    jobname = sprintf('%s_%s', sernum, id);
  end
end
