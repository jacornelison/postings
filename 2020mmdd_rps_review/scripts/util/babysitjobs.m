function farmfiles = babysitjobs(nbase_in,resubmit,minJobs,checkOwners,maxQ,varargin)
% function farmfiles = babysitjobs(nbase_in,resubmit,minJobs,checkOwners,maxQ,varargin)
%
% this file was moved from babsitsims to babysitjobs
% 
% Look for .mat farmfiles that are not deleted and are not currently running or
% queued. If requested, resubmit these jobs.
%
% ex. babysitjobs(1205)
%         or
%     babysitjobs(1205,'resubmit');
%         or
%     babysitjobs(1205,'wait10');
%         or
%     babysitjobs([1205,1206],'wait1');
%         or
%     babysitjobs('final','wait1');
%         or
%     babysitjobs('0750/real*.mat','wait1');
%         or
%     babysitjobs({'final', 'c_t'},'wait1');
%         or
%     babysitjobs({0750,1350},'wait1');
%
%     'wait10' means resubmit failed jobs, then wait 10 minutes, check for failed jobs
%     again, then resubmit, then wait 10 minutes, etc.
%  
%  minJobs: (0) stop babysitting before all jobs have finished
%
%  checkOwners: (0) allow for multiple people using the same farm dir, this make's sure that
%               the queues are checked for all owners of the farm files
%             : 'user1,user2' - check the queues for these specific users
%
%  maxQ: fill the queue to a maximum of maxQ jobs, default is many
%  
%  varargin: these are passed forward to farmit. Allows to change queue, maxtime, mem, etc.
%  e.g. babysitjobs(nbase,resubmit,minJobs,checkOwners,maxQ,'mem',20000)
%       babysitjobs(nbase,resubmit,minJobs,checkOwners,maxQ,'maxtime',600)
%  
%  nbase_in can be numeric, or strings defined from the farmfiles dir.
%  

if ~exist('resubmit','var')
  resubmit=[];
end

if ~exist('minJobs','var') | isempty(minJobs)
  minJobs=0;
end

if ~exist('checkOwners','var') | isempty(checkOwners)
  checkOwners=0;
end

% maxQ is the maximum number of jobs that babysitjobs will submit to the
% queue, per default set it to a rediculus number - mean no restriction.
if ~exist('maxQ','var') | isempty(maxQ)
  maxQ=1000;
end

if ~iscell(nbase_in)
  if isstr(nbase_in)
    nbase_in={nbase_in};
  else
    nbase_in = num2cell(nbase_in);
  end
end

notqueued=1;
Njobs=1;
farmfiles={};

% Keep checking jobs if any have failed or any are running
while 1
  Njobs=0;
  for n=1:numel(nbase_in)
      
    nbase=nbase_in{n};
    
    if isnumeric(nbase) 
      nbase=num2str(nbase, '%.4d');
    end
    % These are the .mat farmfiles for the sim, i.e. not yet successfully
    % finished running.
    if isempty(strfind(nbase, '*')) && isempty(strfind(nbase, '?'))
      farmfiles=listfiles(['farmfiles/',nbase,'/*.mat']);
    else
      farmfiles=listfiles(['farmfiles/',nbase]);
    end
  

    if length(farmfiles)==0 continue; end
    
    % This is how many jobs we are still working on
    Njobs=Njobs+numel(farmfiles);
    
    if isstr(checkOwners)
      owners={checkOwners};
    elseif checkOwners
      % these people own the files + the executor of this script:
      owners = get_owners(farmfiles);
    else
      owners{1} = whoami();
    end    
    
    switch whichhost()
      case 'odyssey'
        queuecmd = {'/usr/bin/squeue --noheader -o "%j" -u '};
      otherwise
        error('unsupported system; do not know how to query the cluster')
    end
 
    % fetch the jobs that are currently running:
    result={};
    for owner = owners
      for qcmd = queuecmd
        % jumb over empty entries:
        if strcmp(owner{:},'') continue; end
        if strcmp(owner{:},' ') continue; end

        cmd=[qcmd{:},owner{:}];
        cresult='';
        [status,cresult]=system(cmd);
        % if the command failed, do it again. It is ok to hang
        % here - no joblist no babysitting...
        while status
          pause(10)
          [status,cresult]=system(cmd);
        end
        if ~isempty(cresult)
          cresult = textscan(cresult, '%s\n');
          result = [result; cresult{1}];
        end
      end
    end
    % These are the corresponding job names for the farmfiles
    jobname = cellfunc(@farmfile2jobname, farmfiles);

    % These are the .mat farmfiles that haven't been deleted but are not
    % currently running or queued.
    notqueued = ~ismember(jobname, result);
    farmfiles = farmfiles(notqueued);
    % this is the amount of jobs in the queue
    queued = numel(jobname)-sum(notqueued);

    fprintf('sim %s: %d jobs have failed/are not queued, %d jobs are queued.\n', ...
        nbase, numel(farmfiles), queued);
    if ~strcmp(resubmit,'resubmit') && isempty(strfind(resubmit,'wait'))
      continue
    end

    if maxQ<1e6
      fprintf('Will refill queue to %d jobs.\n', maxQ);
    else
      disp('Will refill queue.');
    end
    % Resubmit these farmfiles
    for i=1:numel(farmfiles)
      if queued >= maxQ
        break
      end
      % check to make sure the job didn't just complete
      if exist_file(farmfiles{i})
        try
          farmit(farmfiles{i},'resubmit',varargin{:});
          queued = queued+1;
        end
      end
    end
        
  end % for loop
  
  if Njobs<=minJobs
    disp(['babysitting finished with Njobs ',num2str(Njobs),' left over'])
    break
  end

  % wait before checking for jobs again
  if strfind(resubmit,'wait')
    disp(['babysitting until Njobs<=',num2str(minJobs),', now at Njobs=',num2str(length(jobname)),' :'])
    nbase_in
    pause(60*str2num(resubmit(5:end)));
  else
    if ~strcmp(resubmit,'resubmit')
      disp('re-run with with "resubmit" or "waitN" option to resubmit jobs');
    end
    break
  end

end % while loop

return

function owners = get_owners(farmfiles)
  owners={};
  for ii = 1:length(farmfiles)
    owners{ii}=get_file_owner(farmfiles{ii});
  end
  % make sure the executor is included:
  owners{end+1}=whoami();
  owners=unique(owners);
return
