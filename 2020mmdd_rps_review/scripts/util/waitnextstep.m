function waitnextstep(jobname, farmstuff, autoflow, freq)
%waitnextstep(jobname, farmstuff, autoflow, freq)
%
%Queries the job queueing system and can pause execution until a job has
%finished.
%
%    jobname      Job name which is searched for in the queue. (The job name
%                 can actually be any regexp-compatible string applied to the
%                 output of 'squeue -o "%.90j"').
%
%    farmstuff    If false, then returns immediately, otherwise wait. 
%
%    autoflow     If true, then the function will return when the flagged job
%                 has finished. If false, then the function will instead pause
%                 and keyboard input will cause the next statement to execute.
%
%    freq         How often (in min) to check for finished jobs. If you have
%                 many thousands of SLURM jobs, don't set freq < 10 min or we
%                 will overload SLURM SQL. Defaults to 1/6 min (10 sec).
%

if ~exist('farmstuff','var') || isempty(farmstuff)
  farmstuff = true;
end

if ~exist('autoflow','var') || isempty(autoflow)
  autoflow = true;
end

if(~exist('freq','var') || isempty(freq))
  freq=1/6;
end

% control flow to next matlab step, use jobname rlz as flag
if farmstuff
  if autoflow
    disp('pausing until the new data are processed...')
    disp('automatically waiting for jobs to finish')
    flag=1;
    while flag
      pause(60*freq);
      username = whoami();
      [stat_slurm,res_slurm]=unix(['squeue -o "%.90j" -u ' username]);  
      if isempty(regexp(res_slurm,jobname));
        flag=0;
      end
    end
  else
    disp('pause until the new data are processed...')
    disp('hit any key when they are finished... check with sacct --state=r,pd')
    pause
  end
end

return
