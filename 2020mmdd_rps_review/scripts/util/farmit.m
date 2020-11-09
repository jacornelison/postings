% farmit(dirname,cmd)
%
% 0. Running an existing job file
%
% To run a job file at the prompt as if it were on the cluster, simply
% provide the job file path as the only argument to farmit, i.e.
%
%   farmit('farmfiles/jobname.mat');
%
% As a special case, this form takes one optional argument 'interact' which
% avoids resetting the debugger settings (i.e. leaves 'dbstop if error'
% active):
%
%   farmit('farmfiles/jobname.mat', 'interact');
%
%
% 1. Submitting a job to the farm.
%
% When you submit a job using farmit, a temporary file is saved with
% information about the job.  Upon successful completion, the temporary
% file is deleted.  Ordinarily, the filename is generated from the
% current system time and a random string.
%
% Take care: farmit can't control the number of jobs that may try
% to run at one time.  This can cause the load on the disks to become
% huge and slow the system down.
%
% Required arguments:
%   - dirname is the directory in which to save
%     the temporary file.  May be a full filename
%     if desired.
%   - cmd is the command to execute.  This may also
%     be a cell array of commands.
%
% Optional arguments are specified as parameter/value pairs, as in
% farmit(...,'parameter',value).
%   - 'overwrite' - Defaults to false. If true, an existing farm file will
%     be overwritten rather than throwing an error. Another option is
%     'skip' where an existing farm file is left untouched.
%   - 'path' - the Matlab search path to use inside the farm job.
%     By default, this is the same as the path in the workspace
%     from which farmit is called.
%   - 'pwd' - the working directory.  By default, the current
%     directory when farmit is called.
%   - 'queue' - the job queue to submit to.  By default,
%     'serial_requeue,kovac'.
%   - 'jobname' - an identifier for the job.  By default, derived
%     from 'BICEP' and the temporary file name.
%   - 'group' - the job group to submit to. Default derived from 'BICEP' and
%     the user name. LSF only.
%   - 'account' - the job account to submit to. Defaults to 'kovac_lab', with
%     'bicepdata_group' valid as well. SLURM only.
%   - 'license' - the license string to use on SLURM.
%   - 'user' - an e-mail address for reports.  By default, your
%     username on Odyssey.
%   - 'outfile', 'errfile' - files for job output and errors.
%     By default, same as temporary file name but with .out, .err.
%   - 'var' - a variable name (or cell array of names) to pass
%     from the current Matlab workspace into the workspace of the
%     farm job.
%   - 'func' - a cell array of function handles (usually subfunctions
%     defined in the calling code) that should be exposed for use in
%     the farmed job.  Example:
%     ..., 'func', {@subfunc1, @subfunc2}, ...
%   - 'v73' - If true, always save Matfiles as -v7.3.
%     This affects matfiles saved inside the farm job, but not
%     the temporary file.
%   - 'mem' - a number that specifies the memory requirement.
%     This is mandatory for slurm and uses a slightly different cmd.
%     on resubmission 'mem','+100' adds 100 MB to the previous spec
%   - 'maxtime' - the maximum run time of the job in minutes.
%     on resubmission 'maxtime','+60' adds 1hr to the previous spec
%   - 'ntasks' - Number of threads/processes (on a single node) to request from
%     the resource scheduler. This also launches matlab-default rather than
%     matlab to remove RC's addition of -singleCompThread to matlab
%     invocations
%   - 'modules' - A cell array of modules to load in the shell environment (or
%     to unload if the name is prepended by '!').
%   - 'args' - extra stuff to tack on to the batch submission command.
%     Can be a string or a cell array.
%   - 'submit' - 0 to create temporary file but not submit.
%   - 'multiworkspace' - 'true' will submit the commands in the cell
%     cmd to multiple Matlab sessions on a single node.  See examples.
%   - 'nice' - Similar to unix nice, allows for optional deprioritization
%     of the job. See the SLURM sbatch manual for specifics.
%
%   - 'usecompiled' - 1 to use compiled version of code instead of .m file.
%     This also prevents Matlab from being run on the farmed node for this job.
%     0 to initialize Matlab on farmed node and execute .m code
%
% Output arguments:
%   When submitting a job, the name of the temporary file is
%   returned.
%
% Examples:
%   >> farmit('.','reduc_initial(''20100723H01_dk113'');', ...
%             'queue','normal_serial');
%
%   >> farmit('.',{'reduc_initial(''20100723H01_dk113'')',...
%                  'reduc_initial(''20100723H02_dk113'')'}, ...
%             'multiworkspace','true');
%
%
% 2. Resubmitting a job.
%
% You can also resubmit a job file: farmit(fname,'resubmit');
% In this case fname may be a cell array, indicating that the
% temporary files are all to be executed in separate Matlab
% processes on the same node.
%  
% the resubmit option will also accept optional arguments. The
% additional argument(s) will change what has been specified
% during the inital submission. To change for instance the 
% the queue during resubmission do
%   >> farmit(fname,'resubmit','queue','normal_serial')
%
%
% 3. Checking job status.
%
% stat=farmit(fname,'status') will return the status of the job.
% The output may be:
%   [] : file not found (job has finished)
%   'RUN','PEND','SUSP' : as indicated by bjobs
%   'OUTOFTIME' : time limit has been exceeded.
%   'NOJOB' : job not listed in scheduler
%
%
% 4. Killing jobs.
%
% farmit(fname,'kill') will attempt to kill the job.
%

function out1=farmit(fname,cmd,varargin)

% Case 1.  Called with one argument, a file.
% load the file and run the appropriate command.
if nargin==1 || isempty(cmd)

  % special-case handling for 'interact' command -- see case 2 below
  if iscell(fname) && length(fname)==2 && strcmp(fname{2},'interact')
    farmit_interact = true;
    fname = fname{1};
  else
    farmit_interact = false;
  end

  if ~exist(fname,'file')
    error('When farmit is called with one argument, it must be a filename.');
  end

  % In this branch, use variable names that are unlikely to match
  % variable names saved by the user in the farm file.  This will
  % prevent conflicts.
  farmit_fname=fname;
  clear fname
  farmit_S=load(farmit_fname);

  if farmit_interact
    % don't let the debug settings be overridden - use the current environment
    farmit_S.dbstatus = dbstatus('-completenames');
  end
  % Set appropriate settings for farm job environment. Store
  % the old settings in farmit_R so we can restore them later.
  farmit_R=farm_environment_settings(farmit_S);

  if isfield(farmit_S,'ontry')
    farmit_S.ontry=farmit_S.ontry+1;
    save(farmit_S.fname,'-STRUCT','farmit_S');
  end

  if ~isempty(farmit_S.var) && ~isempty(farmit_S.val)
    disp(['FARMIT loading saved variables:']);
    for farmit_i=1:length(farmit_S.var)

      % Some variable names are illegal - they would conflict
      % with local variables in use here!
      if ismember(farmit_S.var{farmit_i}, {'farmit_i','farmit_S','farmit_fname'})
        error(['Illegal variable name ' farmit_S.var{farmit_i} ' conflicts with one already used in farmit.']);
      end

      disp(['  ' farmit_S.var{farmit_i}]);
      eval([farmit_S.var{farmit_i} ' = farmit_S.val{farmit_i};']);
    end
  end

  if ~isempty(farmit_S.func) && ~isempty(farmit_S.func_handle)
    disp(['FARMIT restoring saved function handles:']);
    disp(['  Loading function handles from file a second time after path setup.']);
    farmit_S.func_handle=load(farmit_fname,'func_handle');
    farmit_S.func_handle=farmit_S.func_handle.func_handle;

    for farmit_i=1:length(farmit_S.func)
      % Some function names are illegal - they would conflict
      % with local functions in use here!
      if ismember(farmit_S.func{farmit_i}, {'farmit_i','farmit_S','farmit_fname'})
        error(['Illegal function name ' farmit_S.func{farmit_i} ' conflicts with one already used in farmit.']);
      end

      if strcmpi(func2str(farmit_S.func_handle{farmit_i}),'unknown function')
        error(['Failed to find function corresponding to handle for "' farmit_S.func{farmit_i} '".']);
      end

      disp(['  ' farmit_S.func{farmit_i}]);
      eval([farmit_S.func{farmit_i} ' = farmit_S.func_handle{farmit_i};']);
    end
  end

  disp(['FARMIT running farmed commands:']);
  if isfield(farmit_S,'cmd')
    if ~iscell(farmit_S.cmd)
      farmit_S.cmd={farmit_S.cmd};
    end
    for farmit_i=1:length(farmit_S.cmd)
      disp(['> '  farmit_S.cmd{farmit_i}]);

      if farmit_interact
        eval(farmit_S.cmd{farmit_i});
      else
        try
          eval(farmit_S.cmd{farmit_i});
        catch
          err=lasterror;
          disp(['ERROR in farmit job: ' err.message]);
          % Restore workspace settings before dying
          farm_environment_settings(farmit_R);
          % But we never want dbstop on.
          dbclear('all');
          % And barf.
          rethrow(err);
        end
      end
    end
  end

  % Done, so delete the file.
  delete(farmit_fname);

  % Restore settings to previous values. This prevents annoyance when you execute
  % a farmit job in the interactive environment and find your figure settings, etc.
  % screwed up.
  farm_environment_settings(farmit_R);

  return
end


% Case 2.  Called with the file name of an existing temporary file, and
% an instruction to do something with it such as resubmitting.
if nargin>=2 && (~ischar(fname) || ~exist(fname,'dir')) ...
    && any(ismember(cmd, {'resubmit','reset','retry','status','kill','interact'})) ...
    && exist(fname,'file')

  if ~iscell(fname)
    fname = {fname};
  end
  for i=1:length(fname)
    if ~exist(fname{i},'file') && exist([fname{i} '.mat'],'file')
      fname{i}=[fname{i} '.mat'];
    end
  end

  switch(cmd)
    case 'resubmit',
      for i=1:length(fname)
        S{i}=load(fname{i});
        if exist('varargin','var') | ~isnumeric(S{i}.mem)
          S{i}=farmit_set_args(S{i},varargin);
          T = S{i};
          save(T.fname,'-STRUCT','T');
        end
        if isfield(S{i},'ontry')
          S{i}.ontry=0;
          T = S{i};
          save(T.fname,'-STRUCT','T');
        end
      end
      submit_jobfile(S);
    case 'reset',
      for i=1:length(fname)
        S{i}=load(fname{i});
        if exist('varargin','var') | ~isnumeric(S{i}.mem)
          S{i}=farmit_set_args(S{i},varargin);
          T = S{i};
          save(T.fname,'-STRUCT','T');
        end
        if isfield(S{i},'ontry')
          S{i}.ontry=0;
          T = S{i};
          save(T.fname,'-STRUCT','T');
        end
      end
    case 'retry',
      S=load(fname);
      if isfield(S,'ontry')
        if isfield(S,'retry') && ~isempty(S.retry)
          if S.ontry<S.retry
            submit_jobfile(S);
          end
        end
      end 
    case 'status',
      out1=check_status(fname);
    case 'kill',
      kill_job(fname);
    case 'interact'
      if length(fname) > 1
        error('Interactive job processing is only valid for one job file.');
      end
      % call back into case 1 with special syntax
      fname = {fname{1}, 'interact'};
      farmit(fname);
    otherwise,
      error(['Unknown instruction for processing an existing job file: ' cmd]);
  end
  return
end


% Case 3.  Called with (at least) two arguments, a directory and a command.
% Save a matfile in the directory to encapsulate the job, and submit it.
if nargin>=2
  dirname=fname;
  fname='';
  if isempty(dirname)
    dirname='.';
  end
  if ~exist(dirname,'dir') 
    [dirname fname fext]=fileparts(dirname);
    if isempty(fext)
      fext='.mat';
    end
    fname=[fname fext];
    if ~isempty(dirname) && ~exist(dirname,'dir')
      error('When farmit is called with multiple arguments, first must be a directory or file.');
    end
  end
  oldpwd=pwd;
  cd(dirname);
  dirname=pwd;
  cd(oldpwd);
  fname_random = false;
  if isempty(fname)
    fname=[gen_stamp '.mat'];
    fname_random = true;
  end

  S=farmit_set_props(dirname,fname,cmd,varargin);
  % Note: the code that sets S.val *must* be in the main function!
  % The 'evalin' call will not work in a subfunction!
  val={};
  if isfield(S{1},'var') && ~isempty(S{1}.var)
    for i=1:length(S{1}.var)
      try
        val{i}=evalin('caller',S{1}.var{i});
      catch
        if ~ischar(S{1}.var{i})
          error(['Variable names must be of type char.']);
        else
          error(['Unable to get variable "' S{1}.var{i} '" from parent workspace.']);
        end;
      end
    end
  end
  for i=1:length(S)
    T=S{i};
    T.val=val;
    if exist(T.fname)
      if fname_random
        warning('Auto-generated filename collision. Attempting to recover.')
        T.fname = fullfile(dirname, [gen_stamp '.mat']);
        renamecount = 0;
        while exist(T.fname) && renamecount < 5
          T.fname = fullfile(dirname, [gen_stamp '.mat']);
        end
        if renamecount >= 5
          error('Could not recover from filename collision.')
        end
      else
        % Skip over a farm file without updating if that option has been chosen
        if isfield(T,'overwrite') && ischar(T.overwrite) ...
            && strcmp(T.overwrite,'skip')
          fprintf(1, ['Farmfile %s already exists. Continuing without ' ...
              'overwriting.\n'], T.fname)
          continue;
        end

        % Protect a farm file which already exists if no overwrite behavior is
        % specified or if overwriting has been disabled.
        if ~isfield(T,'overwrite') || T.overwrite==false
          error(['Specified farmfile already exists and will not be ' ...
              'overwritten when overwrite==false.'])
        end
      end
    end
    if exist(T.fname,'file') && T.overwrite
      fprintf(1,['Replacing %s with current options since ' ...
          'overwrite==true.\n'], T.fname);
    end
    % Otherwise write the farm file to disk.
    saveandtest(T.fname,'-struct','T');
    % Report when a file was written without being submitted
    if isfield(T,'submit') && T.submit==false
      fprintf(1, 'Created %s.\n', T.fname)
    end
  end

  if(nargout>0)
    if length(S)>1
      tmpS=cell2mat(S);
      out1={tmpS(:).fname};
    else
      out1=S{1}.fname;
    end
  end
  if ~isfield(S{1},'submit') | S{1}.submit~=0;
    submit_jobfile(S);
  end

  return
end

function S=farmit_set_args(S,p)
  % allow for counting of submits by the babysiting by handing
  % is ..,'count_submits',1,.. as args
  % this is then used to ramp up system requirements when with
  % an extra plus in those args, e.g. 'mem','+1000', adding
  % 1000MB on the second submit and so forth
  if ~isfield(S,'count_submits')
    S.('count_submits')=0;
  else
    S.('count_submits')=S.('count_submits')+1;
  end

  % Fill in values specified by user
  % as parameter-value pairs
  for i=1:2:length(p)
    if isstr(p{i+1}) & p{i+1}(1)=='+' & S.('count_submits')>0
      S.(p{i})=S.(p{i}) + str2num(p{i+1}(2:end));
    else
      S.(p{i})=p{i+1};
    end
  end
return

function S=farmit_set_props(dirname,fname,cmd,p)

  % Defaults, can be overridden with
  % parameter-value pairs
  username=whoami();
  S.user=username;
  [dirtmp ftemp fext]=fileparts(fname);
  S.jobname=['BICEP_' ftemp];
  S.group=['/BICEP_' username];
  S.account='kovac_lab';
  S.license=[];
  outpath=dirname;
  S.outfile=outpath;
  S.errfile=outpath;
  S.farmitpath=fileparts(which('farmit'));
  S.queue='serial_requeue,kovac';
  S.pwd=pwd;
  S.userpath=[];
  S.var={};
  S.val={};
  S.func={};
  S.func_handle={};
  S.args=[];
  S.mem=[];
  S.maxtime=3600;
  S.ntasks=1;
  S.modules={};
  S.submit=true;
  S.overwrite=false;
  S.nice=[];
  S.matlabcmd='matlab';

  S = farmit_set_args(S,p);

  % Only store non-default path information
  fullpath = textscan(path(), '%s', 'delimiter',':;');
  dfltpath = textscan(pathdef(), '%s', 'delimiter',':;');
  cup = ~ismember(fullpath{1},dfltpath{1});
  S.userpath = fullpath{1}(cup);

  if isfield(S,'retry') && ~isempty(S.retry) && isfinite(S.retry) && (S.retry>0)
    S.ontry=0;
  end

  if isfield(S,'var') && ~iscell(S.var)
    S.var={S.var};
  end

  if isfield(S,'func') && ~iscell(S.func)
    S.func={S.func};
  end
  if isfield(S,'func')
    for i=1:length(S.func)
      S.func_handle{i}=S.func{i};
      S.func{i}=func2str(S.func{i});
    end
  end

  S.fname=fullfile(dirname,fname);

  if isfield(S,'multiworkspace') && ~isempty(S.multiworkspace) && strcmpi(S.multiworkspace,'true')
    if ~iscell(cmd)
      cmd={cmd};
    end
    S=repmat({S},length(cmd),1);
    for i=1:length(S)
      S{i}.cmd=cmd{i};
    end
    for i=2:length(S)
      S{i}.fname=fullfile(dirname,[gen_stamp '.mat']);
    end
  else
    S.cmd=cmd;
    S={S};
  end

  for i=1:length(S)
    [dirtmp ftemp fext]=fileparts(S{i}.fname);
    if exist(S{i}.outfile,'dir')
      S{i}.outfile=fullfile(S{i}.outfile,[ftemp '.out']);
    end
    if exist(S{i}.errfile,'dir')
      S{i}.errfile=fullfile(S{i}.errfile,[ftemp '.err']);
    end
  end

  return


function submit_jobfile(S)
  host=whichhost;
  if iscell(S)
    if length(S)>1
      if strcmp(host,'odyssey')
        submit_jobfile_odyssey_multiworkspace(S);
      elseif strcmp(host,'itasca')
        submit_jobfile_itasca_multiworkspace(S);
      else
        for i=1:length(S)
          submit_jobfile_spud(S{i});
        end
      end
      return
    else
      S=S{1};
    end
  end
  if(strcmp(host,'odyssey'))
    if isfield(S,'queue') && ismember(S.queue, {'short_serial', ...
        'normal_serial', 'long_serial', 'unrestricted_serial'})
      submit_jobfile_odyssey(S);
    else
      submit_jobfile_slurm(S);
    end
  elseif strcmp(host,'itasca')
    submit_jobfile_itasca(S);
  else
    submit_jobfile_spud(S);
  end

  return

function submit_jobfile_spud(S)
  disp('submit_jobfile_spud')

  system(['matlab -singleCompThread -nodesktop -nosplash -logfile ' S.outfile ' -r ' ...
          '"try; cd(''' S.pwd '''); farmit(''' S.fname ''');' ...
          'catch; x=lasterror; diary(''' S.errfile ''');' ...
          'disp(x.message); for i=1:numel(x.stack);disp(x.stack(i));end;' ...
          'end; exit;" &']);

  return


function submit_jobfile_odyssey(S)

  disp('Odyssey 1.0 is deprecated.')

  scriptstr=['#! /bin/bash\n'];
  if isfield(S,'user') && ~isempty(S.user)
    scriptstr=[scriptstr '#BSUB -u ' S.user '\n'];
  end
  if isfield(S,'jobname') && ~isempty(S.jobname)
    scriptstr=[scriptstr '#BSUB -J ' S.jobname '\n'];
  end
  if isfield(S,'outfile') && ~isempty(S.outfile)
    scriptstr=[scriptstr '#BSUB -o ' S.outfile '\n'];
  end
  if isfield(S,'errfile') && ~isempty(S.errfile)
    scriptstr=[scriptstr '#BSUB -e ' S.errfile '\n'];
  end
  if ~isfield(S,'queue') || isempty(S.queue)
    S.queue='short_serial';
  end
  scriptstr=[scriptstr '#BSUB -q ' S.queue '\n'];
  if isfield(S,'group') && ~isempty(S.group)
    scriptstr=[scriptstr '#BSUB -g ' S.group '\n'];
  end
  if isfield(S,'mem') && ~isempty(S.mem)
    mem_str=sprintf('-R rusage[mem=%d] -M %d',S.mem,S.mem*1048.576);
    scriptstr=[scriptstr '#BSUB ' mem_str '\n'];
  end
  if isfield(S,'maxtime') && ~isempty(S.maxtime)
    maxtime_str=sprintf('-W %d',S.maxtime);
    scriptstr=[scriptstr '#BSUB ' maxtime_str '\n'];
  end
  if isfield(S,'args') && ~isempty(S.args)
    if ~iscell(S.args)
      S.args={S.args};
    end
    for i=1:length(S.args)
      scriptstr=[scriptstr '#BSUB ' S.args{i} '\n'];
    end
  end
  
  % if using the compiled version of code, the queued command needs to be different
  % than typical, and matlab should not be started.  This if statement takes care of
  % that.  It currently only works when calling farmit from runsim, farm_makepairmaps,
  % and farm_coaddpairmaps, where the farmfile is not a random string, and the
  % variables and values are already saved in that farmfile.  It should be possible to
  % get this section of code to work with arbitrary compiled functions where the
  % farmfile is a random string, but that will take more thought, and we may not want
  % to do that, anyway. 
  if isfield(S,'usecompiled') & S.usecompiled
    % set the mcrCache to the /scratch space of the node, otherwise the homedir
    % will be accessed during runtime
    scriptstr=[scriptstr 'TMPDIR=/scratch/.mcrCache4 \n'];
    scriptstr=[scriptstr 'export MCR_CACHE_ROOT=$TMPDIR \n'];
    % get the directory where the compiled program lies.  This should take care of the
    % case that you start matlab from one directory but your analysis pipeline code is
    % in a different directory, added to your matlab path in your startup file.
    compileddir=fileparts(which('reduc_coaddpairmaps.m'));
    % this system call assumes that 
    system_safe(['/bin/echo -en ''' ...
            scriptstr ...
            compileddir '/' S.cmd S.fname '\n' ...
            ' ''' ...
            ' | bsub']);
  % otherwise submit the jobs as usual to use the matlab function
  else
    system_safe(['/bin/echo -en ''' ...
            scriptstr ...
            'matlab -nodesktop -nosplash -nodisplay -r " cd(''\''''' S.farmitpath '''\''''); farmit(''\''''' ...
            S.fname '''\''''); exit; "\n\n'' ' ...
            ' | bsub']);
  end
  return


function submit_jobfile_slurm(S)

  scriptstr=['#! /bin/bash\n'];
  if isfield(S,'user') && ~isempty(S.user)
    scriptstr=[scriptstr '#SBATCH --mail-user=' S.user '\n'];
  end
  if isfield(S,'jobname') && ~isempty(S.jobname)
    scriptstr=[scriptstr '#SBATCH -J ' S.jobname '\n'];
  end
  if isfield(S,'outfile') && ~isempty(S.outfile)
    scriptstr=[scriptstr '#SBATCH -o ' S.outfile '\n'];
  end
  if isfield(S,'errfile') && ~isempty(S.errfile)
    scriptstr=[scriptstr '#SBATCH -e ' S.errfile '\n'];
  end
  % prevent restarting jobs from overwriting the out and error files:
  scriptstr=[scriptstr '#SBATCH --open-mode=append \n'];
  if ~isfield(S,'queue') || isempty(S.queue)
    S.queue='serial_requeue';
  end
  scriptstr=[scriptstr '#SBATCH -p ' S.queue '\n'];
  %scriptstr=[scriptstr '#SBATCH --acctg-freq=60 \n'];
  if isfield(S,'license') && ~isempty(S.license)
    license_str=sprintf('--license=%s',S.license);
    scriptstr=[scriptstr '#SBATCH ' license_str '\n'];
  end
  if isfield(S,'account') && ~isempty(S.account)
    scriptstr=[scriptstr '#SBATCH --account=' S.account '\n'];
  end
  %% FIXME: This should only be temporary until all code is updated
  %%        to request proper amounts of memory.
  if isempty(S.mem)
    disp('[farmit] ''mem'' must be given on slurm. Assuming 8GB.');
    S.mem = 8000;
  end
  if isfield(S,'mem') && ~isempty(S.mem)
%    mem_str=sprintf('--mem-per-cpu=%d',S.mem);
    mem_str=sprintf('--mem=%d',S.mem);
    scriptstr=[scriptstr '#SBATCH ' mem_str '\n'];
  end
  if isfield(S,'maxtime') && ~isempty(S.maxtime)
    maxtime_str=sprintf(' -t %d',S.maxtime);
    scriptstr=[scriptstr '#SBATCH ' maxtime_str '\n'];
  end
  if isfield(S,'ntasks') && ~isempty(S.ntasks) && S.ntasks > 1
    ntasks_str=sprintf('-n %i',S.ntasks);
    scriptstr=[scriptstr '#SBATCH -N 1\n#SBATCH ' ntasks_str '\n'];
  end
  if isfield(S,'nice') && ~isempty(S.nice)
    nice_str=sprintf('--nice=%i',S.nice);
    scriptstr=[scriptstr '#SBATCH ' nice_str '\n'];
  end

  if isfield(S,'args') && ~isempty(S.args)
    if ~iscell(S.args)
      S.args={S.args};
    end
    for i=1:length(S.args)
      scriptstr=[scriptstr '#SBATCH ' S.args{i} '\n'];
    end
  end
  % propagate total memory ulimit to virtual memory.  Avoids SLURM killing the job when
  % Matlab tries to pre-reserve more memory than is yet needed.
  scriptstr=[scriptstr '\nulimit -v `ulimit -m`;\n'];

  % Setup a trap handler for any errors executing the SLURM script.
  % Useful primarily because SLURM doesn't automatically include host
  % information in the output files like LSF did; having this info when
  % asking for help from rchelp is much appreciated.
  scriptstr = [scriptstr ...
    'function errhandler() {\n' ...
    '  echo -e "\\n\\n$1\\nJob failed at $(date)\\n"\n' ...
    '  exec 1>&2\n' ...
    '  echo -e "\\n\\n$1\\nRun Info:"\n' ...
    '  uname -a\n' ...
    '  echo -e "SLURM Job ID $SLURM_JOB_ID"\n' ...
    '  echo\n' ...
    '  exit 1\n' ...
    '}\n' ...
    '\n' ...
    'trap "errhandler \"Bash caught signal SIGINT.\""  INT\n' ...
    'trap "errhandler \"Bash caught signal SIGABRT.\"" ABRT\n' ...
    'trap "errhandler \"Bash caught signal SIGSEGV.\"" SEGV\n' ...
    'trap "errhandler \"Bash caught signal SIGTERM.\"" TERM\n' ...
    '\n'];

  % Also useful to have the hostname printed into every output file anyway.
  scriptstr = [scriptstr ...
    'echo hostname $(hostname) $(date)\n'...
    '\n'];

  if ~isempty(S.modules)
    for ii=1:length(S.modules)
      module = S.modules{ii};
      if module(1) == '!'
        scriptstr = [scriptstr 'module unload ' module(2:end) '\n'];
      else
        scriptstr = [scriptstr 'module load ' module '\n'];
      end
    end
    scriptstr = [scriptstr '\n'];
  end

  % if using the compiled version of code, the queued command needs to be different
  % than typical, and matlab should not be started.  This if statement takes care of
  % that.  It currently only works when calling farmit from runsim, farm_makepairmaps,
  % and farm_coaddpairmaps, where the farmfile is not a random string, and the
  % variables and values are already saved in that farmfile.  It should be possible to
  % get this section of code to work with arbitrary compiled functions where the
  % farmfile is a random string, but that will take more thought, and we may not want
  % to do that, anyway. 
  if isfield(S,'usecompiled') & S.usecompiled
    % set the mcrCache to the /scratch space of the node, otherwise the homedir
    % will be accessed during runtime
    scriptstr=[scriptstr 'TMPDIR=/scratch/.mcrCache4 \n'];
    scriptstr=[scriptstr 'export MCR_CACHE_ROOT=$TMPDIR \n'];   
    % directory where the compiled program lies.
    compileddir=fileparts(which('reduc_coaddpairmaps.m'));

    scriptstr = [scriptstr ...
      'cd ' S.pwd '\n' ...
      compileddir '/' S.cmd S.fname '\n' ...
      '\n'];

  % otherwise submit the jobs as usual to use the matlab function
  else
    % Choose the specific invocation based on whether we're working in
    % single-threaded mode or are requesting multi-threaded processing.
    if isfield(S,'ntasks') && ~isempty(S.ntasks) && S.ntasks > 1
      threadstr = sprintf('maxNumCompThreads(%i); ', S.ntasks);
    else
      threadstr = '';
    end
    cmd = S.matlabcmd;
    % -nojvm is only usable for R2014a and earlier (before the graphics
    % subsystem rewrite).
    if verLessThan('matlab','8.4')
      % On R2009a, using Java graphics is far more unstable, so disable for
      % batch jobs.
      cmd = [cmd ' -nojvm'];
    end
    cmd = [cmd ' -nodesktop -nosplash -nodisplay'];
    % Add the Matlab invocation, including a try-catch statement to force it
    % to exit on fail and not sit in db mode.
    scriptstr = [scriptstr cmd ' -r "' ...
        threadstr ...
        'cd(''\''''' S.farmitpath '''\''''); ' ...
        'try; ' ...
        '  farmit(''\''''' S.fname '''\''''); '...
        'catch err; ' ...
        '  disp(getReport(err)); ' ...
        '  exit(1); ' ...
        'end; ' ...
        'exit(0); '...
      '"\n\n'];
  end

  % Since Matlab traps any system exceptions thrown (SIGSEGV, SIGINT,
  % etc.), we can simply check the exit status to see whether everything
  % ended correctly. If the $? is non-zero, then dump some extra debugging
  % information about the environment.
  scriptstr = [scriptstr ...
    'if [ $? -ne 0 ]; then\n' ...
    '  errhandler "Matlab ended with non-zero exit code."\n' ...
    'fi\n' ...
    'echo -e "\\n\\nJob completed successfully at $(date)\\n"\n']; 

  try
    % Finally, submit whichever script string has been constructed.
    [r,m] = system_safe(['/bin/echo -en ''' scriptstr ''' | sbatch']);
  catch ex
    % Only propagate the error if it's not a socket timed out message since
    % that doesn't usually signal that the job wasn't submitted to the queue,
    % only that the scheduler didn't respond to sbatch fast enough.
    if isempty(strfind(ex.message, 'Socket timed out'))
      rethrow(ex)
    else
      % Still show the output on stderr. This can only approximatley be right
      % since stdout and stderr both get stored in the same string.
      fprintf(2,strtrim(ex.message))
    end
  end
  % We want messages like "Submitted batch job XXXXXXXX" to be printed to
  % the screen. A regular system() call does this, but due to using
  % system_safe() wrapper, we must do it explicitly.
  %
  % Note that if an exception was caught but we proceeded anyway, m will not
  % be defined.
  if exist('m','var')
    disp(strtrim(m))
  end

return


function submit_jobfile_itasca(S)
  disp('submit_jobfile_itasca')

  scriptstr=['#!/bin/bash -l\n'];
  scriptstr=[scriptstr '#PBS -l walltime=00:10:00,nodes=1:ppn=8 \n'];
  scriptstr=[scriptstr '#PBS -m abe \n'];  

  if isfield(S,'jobname') && ~isempty(S.jobname)
    scriptstr=[scriptstr '#PBS -N ' S.jobname '\n'];
  end
  if isfield(S,'outfile') && ~isempty(S.outfile)
    scriptstr=[scriptstr '#PBS -o ' S.outfile '\n'];
  end
  if isfield(S,'errfile') && ~isempty(S.errfile)
    scriptstr=[scriptstr '#PBS -e ' S.errfile '\n'];
  end
  scriptstr=[scriptstr 'cd ' S.pwd '\n'];
  scriptstr=[scriptstr 'module load matlab/R2011b \n'];

  scriptstr = [scriptstr 'matlab -nodesktop -nosplash -nodisplay -r "maxNumCompThreads(8); farmit(''' S.fname '''); exit;" \n ' ];
  file_1 = fopen('qsubtmp.pbs','w');
  fprintf(file_1,scriptstr);
  fclose(file_1);
  [r,m] = system_safe(['qsub qsubtmp.pbs']);
  disp(strtrim(m))
  return


% Submit a single job, but start several Matlab
% workspaces as background processes.  Kind of
% a hack.
function submit_jobfile_odyssey_multiworkspace(S)
  disp('[farmit] mutliworkspace not implemented on slurm.')
  keyboard

%  if ~iscell(S)
%    S={S};
%  end
%  if length(S)<1
%    return
%  end
%  scriptstr=['#!/bin/sh\n\n'];
%  if isfield(S{1},'user') && ~isempty(S{1}.user)
%    scriptstr=[scriptstr '#BSUB -u ' S{1}.user '\n'];
%  end
%  if isfield(S{1},'jobname') && ~isempty(S{1}.jobname)
%    scriptstr=[scriptstr '#BSUB -J ' S{1}.jobname '\n'];
%  end
%  if isfield(S{1},'outfile') && ~isempty(S{1}.outfile)
%    scriptstr=[scriptstr '#BSUB -o ' S{1}.outfile '\n'];
%  end
%  if isfield(S{1},'errfile') && ~isempty(S{1}.errfile)
%    scriptstr=[scriptstr '#BSUB -e ' S{1}.errfile '\n'];
%  end
%  if ~isfield(S{1},'queue') || isempty(S{1}.queue)
%    S{1}.queue='short_serial';
%  end
%  scriptstr=[scriptstr '#BSUB -q ' S{1}.queue '\n'];
%  if isfield(S{1},'group') && ~isempty(S{1}.group)
%    scriptstr=[scriptstr '#BSUB -g ' S{1}.group '\n'];
%  end
%
%  mlstr='';
%  for i=1:length(S)
%    tmpstr=['matlab -nodesktop -nosplash -r \" ' ...
%            'cd(''' S{i}.farmitpath '''); ' ...
%            'farmit(''' S{i}.fname '''); ' ...
%            'exit; \" '];
%    if i>1
%      tmpstr=[tmpstr ' > ' S{i}.outfile ' 2> ' S{i}.errfile];
%    end
%    mlstr=[mlstr tmpstr ' & '];
%  end
%  system_safe(['/bin/echo -en "' ...
%    scriptstr mlstr '\n\n" ' ...
%    ' | bsub']);
  return

% Submit a single job, but start several Matlab
% workspaces as background processes on itasca
function submit_jobfile_itasca_multiworkspace(S)
  if ~iscell(S)
    S={S};
  end
  if length(S)<1
    return
  end

  scriptstr=['#!/bin/bash -l\n'];
  nNodes=2;
  scriptstr=[scriptstr '#PBS -l walltime=00:10:00,nodes=' num2str(nNodes) ':ppn=8 \n'];
  scriptstr=[scriptstr '#PBS -m abe \n'];  
  
  if isfield(S{1},'jobname') && ~isempty(S{1}.jobname)
    scriptstr=[scriptstr '#PBS -N ' S{1}.jobname '\n'];
  end
  if isfield(S{1},'outfile') && ~isempty(S{1}.outfile)
    scriptstr=[scriptstr '#PBS -o ' S{1}.outfile '\n'];
  end
  if isfield(S{1},'errfile') && ~isempty(S{1}.errfile)
    scriptstr=[scriptstr '#PBS -e ' S{1}.errfile '\n'];
  end
  scriptstr=[scriptstr 'cd ' S{1}.pwd '\n'];
  scriptstr=[scriptstr 'module load matlab/R2011b \n'];
  %find how to use the cores of a node completely; 8 cores to be used here hardcoded for itasca
  coreUseMask = getCoreUseMask(length(S),nNodes*8);
  for i=1:length(S)
    cNode = floor(sum(coreUseMask(1:i-1))*(nNodes*8/sum(coreUseMask)));
    scriptstr = [scriptstr 'pbsdsh -o -n ' num2str(cNode) ' matlab -nodesktop -nosplash -nodisplay -r "dbclear all; maxNumCompThreads(' num2str(coreUseMask(i)) '); farmit(''' S{i}.fname '''); exit;" ' ];
%      scriptstr = [scriptstr '/soft/matlab/R2012a/bin/matlab -nodesktop -nosplash -nodisplay -r "maxNumCompThreads(' num2str(coreUseMask(i)) '); farmit(''' S{i}.fname '''); exit;" ' ];
    if i>1
      scriptstr=[scriptstr ' > ' S{i}.outfile ' 2> ' S{i}.errfile];
    end
    scriptstr = [scriptstr ' & \n'];
  end
  scriptstr = [scriptstr 'wait \n'];
  file_1 = fopen('qsubtmp.pbs','w');
  fprintf(file_1,scriptstr);
  fclose(file_1);
  system_safe(['qsub qsubtmp.pbs']);
return

% splits the number of cores as even as possible among the jobs to be run on these cores
% i.e. getCoreUseMask(3,8) -> mask = [2,3,3]
function mask=getCoreUseMask(nJobs,nCores)
  if (nJobs<nCores)
    mask = ones(nJobs,1)*nCores/nJobs;
    mask = round(mask);
    index = 0;
    while (sum(mask) ~= nCores)
      index = index+1;
      if sum(mask)>nCores 
        mask(index) = mask(index)-1;
      end
      if sum(mask)<nCores 
        mask(index) = mask(index)+1;
      end
      index = mod(index,nJobs);
    end
  else
    mask = ones(nJobs,1);
  end
return

% Set appropriate settings for farm job environment. Store
% the old settings in farmit_R so we can restore them later.
function R=farm_environment_settings(S)
  R = [];

  if isfield(S,'pwd')
    R.pwd=pwd;
    cd(S.pwd);
  end

  % .userpath sets up the execution environment before running the command
  % and stores the current path in .path for restoration later.
  if isfield(S,'userpath')
    R.path = path();
    addpath(S.userpath{:})
  end
  % .path is used to restore the path after execution.
  if isfield(S,'path')
    R.path=path;
    path(S.path);
  end

  if isfield(S,'v73')
    R.v73=com.mathworks.services.Prefs.getStringPref('MatfileSaveFormat');
    if islogical(S.v73) && S.v73
      com.mathworks.services.Prefs.setStringPref('MatfileSaveFormat','v7.3');
    elseif ischar(S.v73)
      com.mathworks.services.Prefs.setStringPref('MatfileSaveFormat',S.v73);
    end
  end
 
  R.DefaultFigureVisible=get(0,'DefaultFigureVisible');
  if isfield(S,'DefaultFigureVisible')
    set(0,'DefaultFigureVisible',S.DefaultFigureVisible);
  else
    set(0,'DefaultFigureVisible','off'); 
  end

  R.DefaultFigureRenderer=get(0,'DefaultFigureRenderer');
  if isfield(S,'DefaultFigureRenderer')
    set(0,'DefaultFigureRenderer',S.DefaultFigureRenderer);
  else
    set(0,'DefaultFigureRenderer','zbuffer');
  end

  R.dbstatus=dbstatus('-completenames');
  if isfield(S,'dbstatus')
    dbstop(S.dbstatus);
  else
    dbclear('all');
  end

  return


function stat=check_status(fname)

  host=whichhost;

  if(strcmp(host,'odyssey'))
    stat=check_status_odyssey(fname);
  elseif(strcmp(host,'itasca'))
    stat=check_status_itasca(fname);
  else
    stat=check_status_spud(fname);
  end

  return


function stat=check_status_spud(fname)

  stat=[];

  return


function stat=check_status_odyssey(fname)
  disp('[farmit] status not implemented on slurm.')
  keyboard

%  if ~iscell(fname)
%    fname={fname};
%    output_as_cell=false;
%  else
%    output_as_cell=true;
%  end
%  stat=cell(size(fname));
%  for i=1:length(fname)
%    if ~exist(fname{i},'file')
%      stat{i}=[];
%      % disp(['No job file ' fname '.']);
%      continue
%    end
%    try
%      S=load(fname{i},'jobname');
%    catch
%      stat{i}=[];
%      % disp(['Job file ' fname ' disappeared.']);
%      continue
%    end
%    B=get_bjobs('-a -J',S.jobname,'-WL');
%    if isempty(B)
%      % disp(['Job "' S.jobname '" is not known.']);
%      stat{i}='NOJOB';
%      continue
%    end
%    stat{i}=strtrim(B.stat{1});
%    timeleft=strtrim(B.time_left);
%    if ~strcmp(timeleft,'') && isempty(strfind(timeleft,'L'))
%      stat{i}='OUTOFTIME';
%    end
%    % disp(['Job "' S.jobname '" has status: ' stat{i}]);
%  end
%  if ~output_as_cell
%    stat=stat{1};
%  end

  return

function stat=check_status_itasca(fname)

  if ~iscell(fname)
    fname={fname};
    output_as_cell=false;
  else
    output_as_cell=true;
  end
  stat=cell(size(fname));
  for i=1:length(fname)
    if ~exist(fname{i},'file')
      stat{i}=[];
      % disp(['No job file ' fname '.']);
      continue
    end
    try
      S=load(fname{i},'jobname');
    catch
      stat{i}=[];
      % disp(['Job file ' fname ' disappeared.']);
      continue
    end
%      B=get_bjobs('-a -J',S.jobname,'-WL');
    B=[]
    if isempty(B)
      % disp(['Job "' S.jobname '" is not known.']);
      stat{i}='NOJOB';
      continue
    end
    stat{i}=strtrim(B.stat{1});
    timeleft=strtrim(B.time_left);
    if ~strcmp(timeleft,'') && isempty(strfind(timeleft,'L'))
      stat{i}='OUTOFTIME';
    end
    % disp(['Job "' S.jobname '" has status: ' stat{i}]);
  end
  if ~output_as_cell
    stat=stat{1};
  end

  return


function kill_job(fname)

  host=whichhost;
  
  if(strcmp(host,'odyssey'))
    kill_job_odyssey(fname);
  elseif(strcmp(host,'itasca'))
    kill_job_itasca(fname);
  else
    kill_job_spud(fname);
  end

  return


function kill_job_spud(fname)

  stat=[];

  return


function kill_job_odyssey(fname)
  disp('[farmit] kill not implemented on slurm.')
  keyboard

%  if ~iscell(fname)
%    fname={fname};
%  end
%  for i=1:length(fname)
%    if ~exist(fname{i},'file')
%      disp(['No job file ' fname '.']);
%      continue
%    end
%    try
%      S=load(fname{i},'jobname');
%    catch
%      disp(['Job file ' fname ' disappeared.']);
%      continue
%    end
%    B=get_bjobs('-a -J',S.jobname);
%    if isempty(B)
%      disp(['Job "' S.jobname '" is not known.']);
%      continue
%    end
%    if ~isempty(B.jobid{1})
%      system_safe(['bkill -s 9 ' B.jobid{1}]);
%    end
%  end

  return

function kill_job_itasca(fname)

  if ~iscell(fname)
    fname={fname};
  end
  for i=1:length(fname)
    if ~exist(fname{i},'file')
      disp(['No job file ' fname '.']);
      continue
    end
    try
      S=load(fname{i},'jobname');
    catch
      disp(['Job file ' fname ' disappeared.']);
      continue
    end
%      B=get_bjobs('-a -J',S.jobname);
    B = []

    if isempty(B)
      disp(['Job "' S.jobname '" is not known.']);
      continue
    end
    if ~isempty(B.jobid{1})
%        system_safe(['bkill -s 9 ' B.jobid{1}]);
    end
  end

  return

