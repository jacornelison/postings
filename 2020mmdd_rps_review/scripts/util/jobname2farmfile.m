function farmfile=jobname2farmfile(jobname,verify)
% farmfile=jobname2farmfile(jobname,verify)
%
% If a job name is babysitjobs-compatible, then the corresponding farmfile
% path is constructed.
%
% INPUTS
%   jobname     A job name string. babysitjobs expects job names to be
%               of the form NNNN_XXXY_ID  or BICEP_ID where NNNN is a 4-digit
%               numeric serial number, XXXY is either a simulation number and
%               type or the string literal 'real', ID is any arbitrary string
%               (usually chosen to identify the operation specifically), and
%               BICEP is a string literal.
%
%   verify      Defaults to false. If true, verifies that the job name
%               adhere's to the expected conventions.
%
% OUTPUTS
%   farmfile    A farm file path string. It will be the standard form
%               'farmfiles/NNNN/XXXY_ID.mat' or 'farmfiles/ID.mat'
%               depending on the input form.
%
% SEE ALSO
%   farmfile2jobname
%   farmfilejobname

  if ~exist('verify','var') || isempty(verify)
    verify = false;
  end

  % For the BICEP_* case, be willing to find it anywhere after a path
  % specifier.
  slashidx = find(jobname == '/', 1, 'last');
  if isempty(slashidx)
    slashidx = 0;
  end
  if strcmp(upper(jobname(slashidx+[1:5])), 'BICEP')
    isbicep = true;

    % Grab the path prefix
    if slashidx > 0
      prefix = jobname(1:slashidx);
    else
      prefix = '';
    end
    % Then the rest of the id comes after the 'BICEP_' which is the 7th
    % character after the prefix ends.
    id = jobname((slashidx+7):end);

  % Otherwise, extract necessary information
  else
    isbicep = false;

    % Do a bit more work to verify the serial number form.
    l = strfind(jobname,'_');
    if (numel(l) < 2 || (l(1) ~= 5 || l(2) ~= 10))
      warning('farmfilejobname:nameFormat', ...
          'No serial number could be identified in jobname `%s`', jobname)

      isbicep = true;
    else
      sernum = jobname(1:4);
      simnum = jobname(6:9);
      id     = jobname(11:end);
    end
  end

  if verify
    is4num = @(s) strcmp(s, sprintf('%04d', str2num(s)));

    % Does the serial number match 4 digit number?
    hassernum = isbicep || is4num(sernum);
    % Does the rest of the serial number follow in the id?
    hassimnum = isbicep || strcmp(simnum, 'real') || is4num(simnum);
    % Is there at least something else to identify the jobs?
    hasid = ~isempty(id);

    if ~hassernum || ~hassimnum || ~hasid
      error('farmfilejobname:nameFormatError', ...
          'jobname `%s` does not conform to expected format', jobname)
    end
  end

  if isbicep
    farmfile = sprintf('farmfiles/%s%s.mat', prefix, id);
  else
    farmfile = sprintf('farmfiles/%s/%s_%s.mat', sernum, simnum, id);
  end
end

